#!/bin/bash
# Submit one runAsymmetry job per input HiForest file.
#
# Usage:
#   bash condor/make_condor.sh OUTPUT_DIR --all [--no-submit|-n]
#   bash condor/make_condor.sh OUTPUT_DIR FILELIST.txt [FILELIST.txt ...] [--no-submit|-n]
#
# OUTPUT_DIR   — absolute EOS/AFS path where output ROOT files are written
# --all        — submit every filelist found in data/txt/
# FILELIST.txt — one or more specific filelists to submit
# --no-submit / -n  — generate submission files without submitting
#
# Mode is auto-detected per filelist from the filename:
#   *_DATA_HP*.txt  → --hard-probes
#   *_DATA_ZB*.txt  → --zero-bias
#   *_MC*.txt       → --monte-carlo
#
# Prerequisites:
#   - Build the project first:  cmake --build build  (from repo root)
#   - All five cone JEC files must be present in data/jec/ before submitting
#   - Set CMSSW_SRC in condor/runtime_wrapper.sh

set -euo pipefail

usage() {
    echo "Usage: $0 OUTPUT_DIR --all [--no-submit|-n]" >&2
    echo "       $0 OUTPUT_DIR FILELIST.txt [FILELIST.txt ...] [--no-submit|-n]" >&2
    exit 1
}

if [[ $# -lt 2 ]]; then usage; fi

OUTPUT_DIR="$1"
shift

USE_ALL=false
NO_SUBMIT=false
EXPLICIT_FILELISTS=()

for arg in "$@"; do
    case "${arg}" in
        --all)          USE_ALL=true ;;
        --no-submit|-n) NO_SUBMIT=true ;;
        *.txt)          EXPLICIT_FILELISTS+=("${arg}") ;;
        *) echo "Unknown argument: ${arg}" >&2; usage ;;
    esac
done

if [[ "${USE_ALL}" == false && ${#EXPLICIT_FILELISTS[@]} -eq 0 ]]; then
    echo "ERROR: specify --all or one or more filelist.txt arguments" >&2
    usage
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDOR_DIR="${REPO_ROOT}/condor"
BINARY="${REPO_ROOT}/bin/runAsymmetry"
LIBRARY="${REPO_ROOT}/lib/libl2residuals.so"
DATA_DIR="${REPO_ROOT}/data"
FILELIST_DIR="${REPO_ROOT}/data/txt"

if [[ ! -f "${BINARY}" ]]; then
    echo "ERROR: ${BINARY} not found — build the project first (cmake --build build)" >&2
    exit 1
fi
if [[ ! -f "${LIBRARY}" ]]; then
    echo "ERROR: ${LIBRARY} not found — build the project first" >&2
    exit 1
fi

if [[ "${USE_ALL}" == true ]]; then
    if [[ ! -d "${FILELIST_DIR}" ]]; then
        echo "ERROR: ${FILELIST_DIR} not found" >&2
        exit 1
    fi
    FILELISTS=("${FILELIST_DIR}"/filelist_HiForest_2024ppref*.txt)
    if [[ ! -f "${FILELISTS[0]}" ]]; then
        echo "ERROR: no filelists found in ${FILELIST_DIR}" >&2
        exit 1
    fi
else
    for f in "${EXPLICIT_FILELISTS[@]}"; do
        if [[ ! -f "${f}" ]]; then
            echo "ERROR: filelist not found: ${f}" >&2
            exit 1
        fi
    done
    FILELISTS=("${EXPLICIT_FILELISTS[@]}")
fi

TODAY=$(date +"%Y-%m-%d_%H-%M-%S")
SUBMISSIONS_DIR="${CONDOR_DIR}/submissions"
mkdir -p "${SUBMISSIONS_DIR}"
WORKDIR="${SUBMISSIONS_DIR}/condor_${TODAY}"
mkdir -p "${WORKDIR}"

(
    cd "${WORKDIR}"

    cp "${CONDOR_DIR}/runtime_wrapper.sh" .
    cp "${BINARY}"  runAsymmetry
    cp "${LIBRARY}" libl2residuals.so
    cp -r "${DATA_DIR}" data
    chmod +x runtime_wrapper.sh runAsymmetry

    mkdir -p "${OUTPUT_DIR}"

    TOTAL_JOBS=0
    TOTAL_LISTS=0

    for FILELIST_PATH in "${FILELISTS[@]}"; do
        BASENAME=$(basename "${FILELIST_PATH}" .txt)

        # Auto-detect mode from filename
        if [[ "${BASENAME}" == *_DATA_HP* ]]; then
            MODE="--hard-probes"
        elif [[ "${BASENAME}" == *_DATA_ZB* ]]; then
            MODE="--zero-bias"
        elif [[ "${BASENAME}" == *_MC* ]]; then
            MODE="--monte-carlo"
        else
            echo "  SKIP  ${BASENAME} — cannot infer mode from filename" >&2
            continue
        fi

        # Label: last _-separated token (HP0, ZB3, MC, …)
        LABEL=$(echo "${BASENAME}" | rev | cut -d_ -f1 | rev)

        mkdir -p "logs/${LABEL}/out" "logs/${LABEL}/err" "logs/${LABEL}/log"
        cp "${FILELIST_PATH}" "filelist_${LABEL}.txt"

        SUBMIT_FILE="submit_${LABEL}.condor"
        COUNT=0

        cat > "${SUBMIT_FILE}" <<EOF
Universe                = vanilla
Executable              = $(pwd)/runtime_wrapper.sh

+JobFlavour             = "longlunch"

should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
Transfer_Output_Files   = ""

Transfer_Input_Files    = $(pwd)/runtime_wrapper.sh,$(pwd)/runAsymmetry,$(pwd)/libl2residuals.so,$(pwd)/data

request_cpus            = 1

EOF

        while IFS= read -r INPUT_FILE; do
            [[ -z "${INPUT_FILE}" ]] && continue

            OUTPUT_FILE="${OUTPUT_DIR}/${LABEL}/output_${COUNT}.root"
            mkdir -p "${OUTPUT_DIR}/${LABEL}"

            cat >> "${SUBMIT_FILE}" <<EOF
Arguments = runAsymmetry ${INPUT_FILE} ${OUTPUT_FILE} ${MODE}
Output    = $(pwd)/logs/${LABEL}/out/job_${COUNT}.out
Error     = $(pwd)/logs/${LABEL}/err/job_${COUNT}.err
Log       = $(pwd)/logs/${LABEL}/log/job_${COUNT}.log
Queue

EOF
            COUNT=$((COUNT + 1))
        done < "filelist_${LABEL}.txt"

        if [[ "${NO_SUBMIT}" == true ]]; then
            echo "  ${LABEL} (${MODE}): ${COUNT} jobs → $(pwd)/${SUBMIT_FILE}"
        else
            echo "  Submitting ${LABEL} (${MODE}): ${COUNT} jobs..."
            condor_submit "${SUBMIT_FILE}"
        fi

        TOTAL_JOBS=$((TOTAL_JOBS + COUNT))
        TOTAL_LISTS=$((TOTAL_LISTS + 1))
    done

    echo ""
    echo "Total: ${TOTAL_JOBS} jobs across ${TOTAL_LISTS} filelists"
    if [[ "${NO_SUBMIT}" == false ]]; then
        echo "Output directory: ${OUTPUT_DIR}"
    fi
    echo "Working directory: ${WORKDIR}"
)

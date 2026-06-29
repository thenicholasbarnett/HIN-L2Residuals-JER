#!/bin/bash
# Submit one runAsymmetry job per input HiForest file across all matching filelists.
#
# Usage:
#   bash condor/make_condor.sh JOBNAME OUTPUT_DIR [MODE] [--no-submit|-n]
#
# JOBNAME    — label used in directory and file names
# OUTPUT_DIR — absolute EOS/AFS path where output ROOT files are written
# MODE       — --hard-probes (default) | --monte-carlo | --zero-bias
# --no-submit / -n  — generate submission files without submitting
#
# All filelists in data/txt/ matching the mode pattern are submitted as
# separate Condor batches (one submit file per filelist). The mode-to-glob
# mapping mirrors cfg/2024ppRef.h:
#   --hard-probes  → filelist_HiForest_2024ppref_DATA_HP*.txt  (~5 datasets)
#   --zero-bias    → filelist_HiForest_2024ppref_DATA_ZB*.txt  (~15 datasets)
#   --monte-carlo  → filelist_HiForest_2024ppref_MC*.txt
#
# Prerequisites:
#   - Build the project first:  mkdir build && cd build && cmake .. && make
#   - All five cone JEC files must be present in data/jec/ before submitting

set -euo pipefail

usage() {
    echo "Usage: $0 JOBNAME OUTPUT_DIR [--monte-carlo|--zero-bias|--hard-probes] [--no-submit|-n]" >&2
    exit 1
}

if [[ $# -lt 2 ]]; then usage; fi

JOBNAME="$1"
OUTPUT_DIR="$2"
MODE="--hard-probes"
NO_SUBMIT=false

# Parse optional arguments
shift 2
for arg in "$@"; do
    case "${arg}" in
        --monte-carlo|--zero-bias|--hard-probes) MODE="${arg}" ;;
        --no-submit|-n)                          NO_SUBMIT=true ;;
        *) echo "Unknown argument: ${arg}" >&2; usage ;;
    esac
done

# Locate the repo root (parent of the condor/ directory this script lives in)
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDOR_DIR="${REPO_ROOT}/condor"
BINARY="${REPO_ROOT}/bin/runAsymmetry"
LIBRARY="${REPO_ROOT}/lib/libl2residuals.so"
DATA_DIR="${REPO_ROOT}/data"
FILELIST_DIR="${REPO_ROOT}/data/txt"

# Glob pattern per mode — mirrors cfg/2024ppRef.h naming convention
case "${MODE}" in
    --hard-probes)  FILELIST_GLOB="filelist_HiForest_2024ppref_DATA_HP*.txt" ;;
    --zero-bias)    FILELIST_GLOB="filelist_HiForest_2024ppref_DATA_ZB*.txt" ;;
    --monte-carlo)  FILELIST_GLOB="filelist_HiForest_2024ppref_MC*.txt" ;;
esac

# Sanity checks
if [[ ! -f "${BINARY}" ]]; then
    echo "ERROR: ${BINARY} not found — build the project first (cmake .. && make)" >&2
    exit 1
fi
if [[ ! -f "${LIBRARY}" ]]; then
    echo "ERROR: ${LIBRARY} not found — build the project first" >&2
    exit 1
fi
if [[ ! -d "${DATA_DIR}" ]]; then
    echo "ERROR: ${DATA_DIR} not found" >&2
    exit 1
fi

# Collect matching filelists
FILELISTS=("${FILELIST_DIR}"/${FILELIST_GLOB})
if [[ ! -f "${FILELISTS[0]}" ]]; then
    echo "ERROR: no filelists matched ${FILELIST_DIR}/${FILELIST_GLOB}" >&2
    exit 1
fi

TODAY=$(date +"%Y-%m-%d_%H-%M-%S")
WORKDIR="$(pwd)/condor_${JOBNAME}_${TODAY}"
mkdir -p "${WORKDIR}"

(
    cd "${WORKDIR}"

    # Copy shared artifacts once — all submit files in this work dir reference them.
    cp "${CONDOR_DIR}/runtime_wrapper.sh" .
    cp "${BINARY}"  runAsymmetry
    cp "${LIBRARY}" libl2residuals.so
    cp -r "${DATA_DIR}" data
    chmod +x runtime_wrapper.sh runAsymmetry

    mkdir -p "${OUTPUT_DIR}"

    TOTAL_JOBS=0

    for FILELIST_PATH in "${FILELISTS[@]}"; do
        # Extract dataset label: HP0, ZB3, MC, etc. (last _-separated token, no extension)
        LABEL=$(basename "${FILELIST_PATH}" .txt | rev | cut -d_ -f1 | rev)

        mkdir -p "logs/${LABEL}/out" "logs/${LABEL}/err" "logs/${LABEL}/log"
        cp "${FILELIST_PATH}" "filelist_${LABEL}.txt"

        SUBMIT_FILE="submit_${JOBNAME}_${LABEL}.condor"
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
            echo "  ${LABEL}: ${COUNT} jobs written to $(pwd)/${SUBMIT_FILE}"
        else
            echo "  Submitting ${LABEL}: ${COUNT} jobs..."
            condor_submit "${SUBMIT_FILE}"
        fi

        TOTAL_JOBS=$((TOTAL_JOBS + COUNT))
    done

    echo ""
    echo "Total jobs: ${TOTAL_JOBS} across ${#FILELISTS[@]} filelists (mode: ${MODE})"
    if [[ "${NO_SUBMIT}" == false ]]; then
        echo "Output directory: ${OUTPUT_DIR}"
    fi
    echo "Working directory: ${WORKDIR}"
)

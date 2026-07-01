#!/bin/bash
# Submit one runAsymmetry job per input HiForest file.
#
# Usage:
#   bash condor/make_condor.sh OUTPUT_DIR -a|--all-txt [--no-submit|-n] [TAG=value] [CONFIG=path]
#   bash condor/make_condor.sh OUTPUT_DIR FILELIST.txt [FILELIST.txt ...] [--no-submit|-n] [TAG=value] [CONFIG=path]
#
# OUTPUT_DIR    — absolute EOS/AFS path where output ROOT files are written
# -a/--all-txt  — submit every .txt filelist found in data/txt/
# FILELIST.txt  — one or more specific filelists to submit
# --no-submit/-n — generate submission files without submitting
# TAG=value     — optional label for this pass (e.g. TAG=abs_eta, TAG=clos_dir_eta).
#                 Output goes to OUTPUT_DIR/condor/asymmetry_<value>/<timestamp>
#                 instead of OUTPUT_DIR/condor/asymmetry/<timestamp> — use this to
#                 keep separate closure/reprocessing passes from landing in the
#                 same output tree. No slashes allowed in value. Independent of
#                 CONFIG below — mix and match freely, neither implies the other.
# CONFIG=path   — TOML to submit with (default: cfg/2024ppRef.toml). Whatever is
#                 selected is transferred to the sandbox under a fixed name
#                 (analysis_config.toml), so runtime_wrapper.sh never needs to
#                 know the source filename — pointing this at a different run
#                 period/collision system's TOML needs no other change here.
#
# Mode is auto-detected per filelist from the filename (case-insensitive):
#   *hp* or *hardprobes*   → --hard-probes
#   *zb* or *zerobias*     → --zero-bias
#   *mc* or *montecarlo*   → --monte-carlo
#
# Prerequisites:
#   - Build the project first:  cmake --build build  (from repo root)
#   - All five cone JEC files must be present in data/jec/ before submitting
#   - Set CMSSW_SRC in condor/runtime_wrapper.sh

set -euo pipefail

usage() {
    echo "Usage: $0 OUTPUT_DIR -a|--all-txt [--no-submit|-n] [TAG=value] [CONFIG=path]" >&2
    echo "       $0 OUTPUT_DIR FILELIST.txt [FILELIST.txt ...] [--no-submit|-n] [TAG=value] [CONFIG=path]" >&2
    exit 1
}

if [[ $# -lt 2 ]]; then usage; fi

OUTPUT_DIR="$1"
shift

USE_ALL=false
NO_SUBMIT=false
EXPLICIT_FILELISTS=()
RUN_TAG=""
CONFIG_PATH=""

for arg in "$@"; do
    case "${arg}" in
        -a|--all-txt)   USE_ALL=true ;;
        --no-submit|-n) NO_SUBMIT=true ;;
        TAG=*)          RUN_TAG="${arg#TAG=}" ;;
        CONFIG=*)       CONFIG_PATH="${arg#CONFIG=}" ;;
        *.txt)          EXPLICIT_FILELISTS+=("${arg}") ;;
        *) echo "Unknown argument: ${arg}" >&2; usage ;;
    esac
done

if [[ "${RUN_TAG}" == */* ]]; then
    echo "ERROR: TAG must not contain '/': ${RUN_TAG}" >&2
    exit 1
fi

if [[ "${USE_ALL}" == false && ${#EXPLICIT_FILELISTS[@]} -eq 0 ]]; then
    echo "ERROR: specify -a/--all-txt or one or more filelist.txt arguments" >&2
    usage
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDOR_DIR="${REPO_ROOT}/condor"
BINARY="${REPO_ROOT}/bin/runAsymmetry"
LIBRARY="${REPO_ROOT}/lib/libl2residuals.so"
DATA_DIR="${REPO_ROOT}/data"
FILELIST_DIR="${REPO_ROOT}/data/txt"

if [[ -z "${CONFIG_PATH}" ]]; then CONFIG_PATH="${REPO_ROOT}/cfg/2024ppRef.toml"; fi
if [[ ! -f "${CONFIG_PATH}" ]]; then
    echo "ERROR: CONFIG file not found: ${CONFIG_PATH}" >&2
    exit 1
fi
CONFIG_PATH="$(cd "$(dirname "${CONFIG_PATH}")" && pwd)/$(basename "${CONFIG_PATH}")"

if [[ -z "${CMSSW_BASE:-}" ]]; then
    echo "ERROR: cmsenv is not active. The binary must be built and submitted from a cmsenv shell." >&2
    echo "       source /cvmfs/cms.cern.ch/cmsset_default.sh" >&2
    echo "       cd <CMSSW>/src && cmsenv && cd -" >&2
    echo "       cmake --build build  # if not already built" >&2
    exit 1
fi

if [[ ! -f "${BINARY}" ]]; then
    echo "ERROR: ${BINARY} not found — run: cmake --build build" >&2
    exit 1
fi
if [[ ! -f "${LIBRARY}" ]]; then
    echo "ERROR: ${LIBRARY} not found — run: cmake --build build" >&2
    exit 1
fi

# Belt-and-suspenders: verify the binary's embedded RPATH points to cvmfs.
# Catches the case where cmsenv was sourced after the binary was built without it.
BINARY_RPATH="$(readelf -d "${BINARY}" 2>/dev/null | grep -E 'RPATH|RUNPATH' | grep -o '\[.*\]' | tr -d '[]' || true)"
if [[ -n "${BINARY_RPATH}" ]] && ! echo "${BINARY_RPATH}" | grep -q '/cvmfs/'; then
    echo "ERROR: ${BINARY} RPATH does not point to cvmfs: ${BINARY_RPATH}" >&2
    echo "       Binary was built without cmsenv. Rebuild: cmake --build build" >&2
    exit 1
fi

if [[ "${USE_ALL}" == true ]]; then
    if [[ ! -d "${FILELIST_DIR}" ]]; then
        echo "ERROR: ${FILELIST_DIR} not found" >&2
        exit 1
    fi
    FILELISTS=("${FILELIST_DIR}"/*.txt)
    if [[ ! -f "${FILELISTS[0]}" ]]; then
        echo "ERROR: no .txt filelists found in ${FILELIST_DIR}" >&2
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
WORKDIR="${SUBMISSIONS_DIR}/${TODAY}"
mkdir -p "${WORKDIR}"
# asymmetry_<TAG> keeps separate passes (e.g. TAG=abs_eta vs TAG=clos_dir_eta)
# from landing in the same output tree; plain "asymmetry" when no TAG is given.
ASYM_DIR_NAME="asymmetry"
if [[ -n "${RUN_TAG}" ]]; then ASYM_DIR_NAME="asymmetry_${RUN_TAG}"; fi

# Normalize: strip trailing /condor/<asym_dir_name> or /condor so the user can pass
# any of outdir, outdir/condor, or outdir/condor/<asym_dir_name> and land in the same place.
OUTPUT_DIR="${OUTPUT_DIR%/}"
OUTPUT_DIR="${OUTPUT_DIR%/condor/${ASYM_DIR_NAME}}"
OUTPUT_DIR="${OUTPUT_DIR%/condor}"
OUTPUT_DIR="${OUTPUT_DIR}/condor/${ASYM_DIR_NAME}/${TODAY}"

source "$(dirname "${BASH_SOURCE[0]}")/draw_bar.sh"

(
    cd "${WORKDIR}"

    cp "${CONDOR_DIR}/runtime_wrapper.sh" .
    cp "${BINARY}"  runAsymmetry
    cp "${LIBRARY}" libl2residuals.so
    cp -r "${DATA_DIR}" data
    cp "${CONFIG_PATH}" analysis_config.toml
    chmod +x runtime_wrapper.sh runAsymmetry

    mkdir -p "${OUTPUT_DIR}"

    TOTAL_JOBS=0
    TOTAL_LISTS=0

    for FILELIST_PATH in "${FILELISTS[@]}"; do
        BASENAME=$(basename "${FILELIST_PATH}" .txt)

        # Auto-detect mode from filename (case-insensitive; matches hp/hardprobes, zb/zerobias, mc/montecarlo)
        LOWER="${BASENAME,,}"
        if [[ "${LOWER}" == *hp* || "${LOWER}" == *hardprobes* ]]; then
            MODE="--hard-probes"
        elif [[ "${LOWER}" == *zb* || "${LOWER}" == *zerobias* ]]; then
            MODE="--zero-bias"
        elif [[ "${LOWER}" == *mc* || "${LOWER}" == *montecarlo* ]]; then
            MODE="--monte-carlo"
        else
            echo "  SKIP  ${BASENAME} — cannot infer mode from filename" >&2
            continue
        fi

        # Label: last _-separated token (HP0, ZB3, MC, …)
        LABEL=$(echo "${BASENAME}" | rev | cut -d_ -f1 | rev)

        mkdir -p "logs/${LABEL}/out" "logs/${LABEL}/err" "logs/${LABEL}/log"

        FILELIST_FILE="data/txt/$(basename "${FILELIST_PATH}")"
        TOTAL=$(grep -c . "${FILELIST_FILE}" || echo 0)
        SUBMIT_FILE="submit_${LABEL}.condor"
        COUNT=0

        cat > "${SUBMIT_FILE}" <<EOF
Universe                = vanilla
Executable              = $(pwd)/runtime_wrapper.sh

+JobFlavour             = "longlunch"
JobBatchName            = "${LABEL}"

should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
Transfer_Output_Files   = ""

Transfer_Input_Files    = $(pwd)/runtime_wrapper.sh,$(pwd)/runAsymmetry,$(pwd)/libl2residuals.so,$(pwd)/data,$(pwd)/analysis_config.toml

request_cpus            = 1

EOF

        pick_bar_color
        start_bar_timer
        draw_bar "${BAR_COLOR}" "${LABEL}:" 0 "${TOTAL}"

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
            draw_bar "${BAR_COLOR}" "${LABEL}:" "${COUNT}" "${TOTAL}"
        done < "${FILELIST_FILE}"

        printf "\n\n"

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
    echo "Config used: ${CONFIG_PATH} (transferred as analysis_config.toml, also archived in ${WORKDIR})"
    if [[ "${NO_SUBMIT}" == false ]]; then
        echo "Output directory: ${OUTPUT_DIR}"
    fi
    echo "Working directory: ${WORKDIR}"
)

#!/bin/bash
# Submit one runAsymmetry job per input HiForest file.
#
# Usage:
#   bash condor/make_condor.sh JOBNAME OUTPUT_DIR [MODE] [--no-submit|-n]
#
# JOBNAME    — label used in directory and file names
# OUTPUT_DIR — absolute EOS/AFS path where output ROOT files are written
# MODE       — --hard-probes (default) | --monte-carlo | --zero-bias
# --no-submit / -n  — generate the submission file without submitting
#
# The input filelist is selected automatically from MODE, mirroring the
# kFilelistHP / kFilelistZB / kFilelistMC constants in cfg/2024ppRef.h:
#   --hard-probes  → data/txt/filelist_HiForest_2024ppref_DATA_HP0.txt
#   --zero-bias    → data/txt/filelist_HiForest_2024ppref_DATA_ZB0.txt
#   --monte-carlo  → data/txt/filelist_HiForest_2024ppref_MC.txt
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

# Mirror kFilelistHP / kFilelistZB / kFilelistMC from cfg/2024ppRef.h
FILELIST_HP="${REPO_ROOT}/data/txt/filelist_HiForest_2024ppref_DATA_HP0.txt"
FILELIST_ZB="${REPO_ROOT}/data/txt/filelist_HiForest_2024ppref_DATA_ZB0.txt"
FILELIST_MC="${REPO_ROOT}/data/txt/filelist_HiForest_2024ppref_MC.txt"

case "${MODE}" in
    --hard-probes)  FILELIST="${FILELIST_HP}" ;;
    --zero-bias)    FILELIST="${FILELIST_ZB}" ;;
    --monte-carlo)  FILELIST="${FILELIST_MC}" ;;
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
if [[ ! -f "${FILELIST}" ]]; then
    echo "ERROR: filelist not found: ${FILELIST}" >&2
    exit 1
fi

TODAY=$(date +"%Y-%m-%d_%H-%M-%S")
WORKDIR="$(pwd)/condor_${JOBNAME}_${TODAY}"
mkdir -p "${WORKDIR}"

(
    cd "${WORKDIR}"
    mkdir -p logs/out logs/err logs/log

    # Copy everything the worker nodes will need into the work directory.
    # Condor transfers these to the sandbox on each worker.
    cp "${CONDOR_DIR}/runtime_wrapper.sh" .
    cp "${BINARY}"  runAsymmetry
    cp "${LIBRARY}" libl2residuals.so
    cp -r "${DATA_DIR}" data
    cp "${FILELIST}" filelist.txt
    chmod +x runtime_wrapper.sh runAsymmetry

    mkdir -p "${OUTPUT_DIR}"

    SUBMIT_FILE="submit_${JOBNAME}.condor"

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

    COUNT=0

    while IFS= read -r INPUT_FILE; do
        [[ -z "${INPUT_FILE}" ]] && continue

        OUTPUT_FILE="${OUTPUT_DIR}/output_${JOBNAME}_${COUNT}.root"

        cat >> "${SUBMIT_FILE}" <<EOF
Arguments = runAsymmetry ${INPUT_FILE} ${OUTPUT_FILE} ${MODE}
Output    = $(pwd)/logs/out/job_${COUNT}.out
Error     = $(pwd)/logs/err/job_${COUNT}.err
Log       = $(pwd)/logs/log/job_${COUNT}.log
Queue

EOF
        COUNT=$((COUNT + 1))

    done < filelist.txt

    if [[ "${NO_SUBMIT}" == true ]]; then
        echo "--no-submit flag set — skipping condor_submit."
        echo "${COUNT} jobs written to $(pwd)/${SUBMIT_FILE}"
    else
        echo "Submitting ${COUNT} jobs (mode: ${MODE})..."
        condor_submit "${SUBMIT_FILE}"
        echo "Done. Output directory: ${OUTPUT_DIR}"
    fi

    echo "Working directory: ${WORKDIR}"
)

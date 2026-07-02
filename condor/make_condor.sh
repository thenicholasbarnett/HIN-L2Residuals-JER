#!/bin/bash
# Submit one runAsymmetry job per input HiForest file.
#
# Usage:
#   bash condor/make_condor.sh OUTPUT_DIR -a|--all-txt CONFIG=path [--no-submit|-n] [TAG=value]
#   bash condor/make_condor.sh OUTPUT_DIR FILELIST.txt [FILELIST.txt ...] CONFIG=path [--no-submit|-n] [TAG=value]
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
# CONFIG=path   — required; TOML to submit with. Which TOML gets submitted is a
#                 physics-affecting choice (e.g. main run-period config vs. a
#                 closure config with different residual_files), so there is no
#                 implicit default — pass CONFIG=cfg/2024ppRef.toml explicitly
#                 for the standard run-period config. Whatever is selected is
#                 transferred to the sandbox under a fixed name
#                 (analysis_config.toml), so runtime_wrapper.sh never needs to
#                 know the source filename — pointing this at a different run
#                 period/collision system's TOML needs no other change here.
#
# Mode for each filelist is looked up from [condor.filelist_modes] in the
# selected TOML, keyed by the filelist's basename (without .txt):
#   [condor.filelist_modes]
#   filelist_HiForest_2024ppref_DATA_HP0 = "triggered"
#   filelist_HiForest_2024ppref_DATA_ZB0 = "non-triggered"
#   filelist_HiForest_2024ppref_MC       = "mc"
# Comment out (or simply omit) a filelist's entry to skip submitting jobs for
# it without moving or renaming the file. Values: triggered | non-triggered | mc.
# The physical dataset name (HP0, ZB0, MinBias, Jet80, ...) is free-form —
# only the mapped mode value matters; the key just has to match the filelist's
# basename exactly.
#
# Prerequisites:
#   - On lxplus, run cmsenv before configuring/building this repo:
#       source /cvmfs/cms.cern.ch/cmsset_default.sh
#       cd <CMSSW_RELEASE>/src && cmsenv && cd /path/to/L2Residuals-2024ppref
#       cmake -S . -B build && cmake --build build
#     or, if the repo is checked out as $CMSSW_BASE/src/L2Residuals:
#       scram b -j4
#   - All five cone JEC files must be present in data/jec/ before submitting
#   - Set [condor].cmssw_src in the selected TOML

set -euo pipefail

usage() {
    echo "Usage: $0 OUTPUT_DIR -a|--all-txt CONFIG=path [--no-submit|-n] [TAG=value]" >&2
    echo "       $0 OUTPUT_DIR FILELIST.txt [FILELIST.txt ...] CONFIG=path [--no-submit|-n] [TAG=value]" >&2
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
CMAKE_BINARY="${REPO_ROOT}/build/bin/runAsymmetry"
CMAKE_LIBRARY="${REPO_ROOT}/build/lib/libl2residuals.so"
SCRAM_BINARY=""
BINARY="${CMAKE_BINARY}"
LIBRARY="${CMAKE_LIBRARY}"
if [[ -n "${CMSSW_BASE:-}" && -n "${SCRAM_ARCH:-}" && -x "${CMSSW_BASE}/bin/${SCRAM_ARCH}/runAsymmetry" ]]; then
    SCRAM_BINARY="${CMSSW_BASE}/bin/${SCRAM_ARCH}/runAsymmetry"
    BINARY="${SCRAM_BINARY}"
    LIBRARY=""
fi
DATA_DIR="${REPO_ROOT}/data"
FILELIST_DIR="${REPO_ROOT}/data/txt"

if [[ -z "${CONFIG_PATH}" ]]; then
    echo "ERROR: CONFIG=path is required -- which TOML gets submitted is a physics-affecting" >&2
    echo "       choice (e.g. main run-period config vs. a closure config), so there is no" >&2
    echo "       implicit default. Pass CONFIG=cfg/2024ppRef.toml explicitly if that's what you want." >&2
    exit 1
fi
if [[ ! -f "${CONFIG_PATH}" ]]; then
    echo "ERROR: CONFIG file not found: ${CONFIG_PATH}" >&2
    exit 1
fi
CONFIG_PATH="$(cd "$(dirname "${CONFIG_PATH}")" && pwd)/$(basename "${CONFIG_PATH}")"

# Look up the run mode for a filelist basename from [condor.filelist_modes] in
# the selected TOML. Prints "triggered", "non-triggered", "mc", or nothing if
# there's no active (uncommented) entry for that key.
lookup_filelist_mode() {
    local key="$1"
    awk -v key="${key}" '
        /^[[:space:]]*\[/ {
            insec = ( $0 ~ /^[[:space:]]*\[condor\.filelist_modes\][[:space:]]*(#.*)?$/ )
            next
        }
        insec {
            line = $0
            sub( /#.*/, "", line )
            if ( line ~ ( "^[[:space:]]*" key "[[:space:]]*=" ) ) {
                val = line
                sub( /^[^=]*=[[:space:]]*"/, "", val )
                sub( /".*/, "", val )
                print val
                exit
            }
        }
    ' "${CONFIG_PATH}"
}

CMSSW_SRC_FROM_CONFIG="$(awk '
    /^[[:space:]]*\[/ {
        in_condor = ( $0 ~ /^[[:space:]]*\[condor\][[:space:]]*(#.*)?$/ )
        next
    }
    in_condor && $0 ~ /^[[:space:]]*cmssw_src[[:space:]]*=/ {
        line = $0
        sub( /^[^=]*=[[:space:]]*"/, "", line )
        sub( /"[[:space:]]*(#.*)?$/, "", line )
        print line
        exit
    }
' "${CONFIG_PATH}")"
if [[ -z "${CMSSW_SRC_FROM_CONFIG}" || "${CMSSW_SRC_FROM_CONFIG}" == REQUIRED_* ]]; then
    echo "ERROR: ${CONFIG_PATH} must set [condor].cmssw_src to the CMSSW src directory used for worker jobs" >&2
    exit 1
fi
if [[ ! -d "${CMSSW_SRC_FROM_CONFIG}" ]]; then
    echo "ERROR: [condor].cmssw_src does not exist: ${CMSSW_SRC_FROM_CONFIG}" >&2
    exit 1
fi

if [[ -z "${CMSSW_BASE:-}" ]]; then
    echo "ERROR: cmsenv is not active. The binary must be built and submitted from a cmsenv shell." >&2
    echo "       source /cvmfs/cms.cern.ch/cmsset_default.sh" >&2
    echo "       cd <CMSSW>/src && cmsenv && cd -" >&2
    echo "       cd /path/to/L2Residuals-2024ppref" >&2
    echo "       cmake -S . -B build && cmake --build build" >&2
    echo "       or build the SCRAM executable with: scram b -j4" >&2
    exit 1
fi

if [[ ! -f "${BINARY}" ]]; then
    echo "ERROR: no runAsymmetry executable found." >&2
    echo "       For CMake, run: cmake --build build" >&2
    echo "       For SCRAM, put the repo at $CMSSW_BASE/src/L2Residuals and run: scram b -j4" >&2
    exit 1
fi
if [[ -n "${LIBRARY}" && ! -f "${LIBRARY}" ]]; then
    echo "ERROR: ${LIBRARY} not found — run: cmake --build build" >&2
    exit 1
fi

# Belt-and-suspenders: verify the binary's embedded RPATH points to cvmfs.
# Catches the case where cmsenv was sourced after the binary was built without it.
BINARY_RPATH="$(readelf -d "${BINARY}" 2>/dev/null | grep -E 'RPATH|RUNPATH' | grep -o '\[.*\]' | tr -d '[]' || true)"
if [[ -z "${SCRAM_BINARY}" && -n "${BINARY_RPATH}" ]] && ! echo "${BINARY_RPATH}" | grep -q '/cvmfs/'; then
    echo "ERROR: ${BINARY} RPATH does not point to cvmfs: ${BINARY_RPATH}" >&2
    echo "       Binary was built without cmsenv. Reconfigure and rebuild after cmsenv:" >&2
    echo "       cmake -S . -B build && cmake --build build" >&2
    echo "       or build the SCRAM executable with: scram b -j4" >&2
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

    CMSSW_SRC_ESCAPED="${CMSSW_SRC_FROM_CONFIG//\\/\\\\}"
    CMSSW_SRC_ESCAPED="${CMSSW_SRC_ESCAPED//&/\\&}"
    sed "s|@CMSSW_SRC@|${CMSSW_SRC_ESCAPED}|g" "${CONDOR_DIR}/runtime_wrapper.sh" > runtime_wrapper.sh
    cp "${BINARY}"  runAsymmetry
    if [[ -n "${LIBRARY}" ]]; then cp "${LIBRARY}" libl2residuals.so; fi
    cp -r "${DATA_DIR}" data
    cp "${CONFIG_PATH}" analysis_config.toml
    chmod +x runtime_wrapper.sh runAsymmetry

    mkdir -p "${OUTPUT_DIR}"

    TOTAL_JOBS=0
    TOTAL_LISTS=0

    for FILELIST_PATH in "${FILELISTS[@]}"; do
        BASENAME=$(basename "${FILELIST_PATH}" .txt)

        FILELIST_MODE="$(lookup_filelist_mode "${BASENAME}")"
        case "${FILELIST_MODE}" in
            triggered|non-triggered|mc)
                MODE="${FILELIST_MODE}"
                ;;
            "")
                echo "  SKIP  ${BASENAME} — no [condor.filelist_modes] entry in ${CONFIG_PATH} (missing or commented out)" >&2
                continue
                ;;
            *)
                echo "  SKIP  ${BASENAME} — unknown mode \"${FILELIST_MODE}\" in [condor.filelist_modes] (expected triggered|non-triggered|mc)" >&2
                continue
                ;;
        esac

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

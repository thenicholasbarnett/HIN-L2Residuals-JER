#!/bin/bash
# Submit one runAsymmetry job per input HiForest file.
#
# Arguments are JetMET-style "-flag value" (matching the compiled binaries'
# CLI, external/jetmet/CommandLine.h) or bare boolean "-flag" switches -- no
# positional arguments, no KEY=value. Flag names are case-insensitive
# (matched after uppercasing) but otherwise this is a plain left-to-right
# parse: an unrecognized flag, a value-taking flag with nothing after it, or
# any argument that isn't a "-flag" at all is an immediate usage error.
#
# Usage:
#   bash condor/make_condor.sh -output dir -alltxt -config path [-nosubmit] [-tag value]
#   bash condor/make_condor.sh -output dir -filelists a.txt b.txt ... -config path [-nosubmit] [-tag value]
#
# -output dir       — required; absolute EOS/AFS path where output ROOT files are written
# -alltxt           — bare switch; submit every .txt filelist found in data/txt/ (default off)
# -filelists ...    — one or more filelists and/or directories to submit, as separate
#                      words (e.g. -filelists a.txt b.txt, or -filelists data/txt/2026-07-03).
#                      Consumes every following word up to the next "-flag". A directory
#                      entry expands to every *.txt file directly inside it (non-recursive),
#                      so you can stage a subset of filelists into its own directory (e.g.
#                      named after today's date) and point -filelists at that directory
#                      instead of listing each file. [condor.filelist_modes] keys off the
#                      filelist's basename only, not its path, so this works with no TOML
#                      changes regardless of which directory a filelist lives in. Mutually
#                      exclusive with -alltxt; exactly one of the two is required.
# -nosubmit         — bare switch; generate submission files without submitting (default off)
# -tag value        — optional label for this pass (e.g. -tag abs_eta, -tag clos_dir_eta).
#                      Output goes to OUTPUT_DIR/condor/asymmetry_<value>/<timestamp>
#                      instead of OUTPUT_DIR/condor/asymmetry/<timestamp> — use this to
#                      keep separate closure/reprocessing passes from landing in the
#                      same output tree. No slashes allowed in value. Independent of
#                      -config below — mix and match freely, neither implies the other.
# -config path      — required; TOML to submit with. Which TOML gets submitted is a
#                      physics-affecting choice (e.g. main run-period config vs. a
#                      closure config with different residual_files), so there is no
#                      implicit default — pass -config cfg/2024ppRef.toml explicitly
#                      for the standard run-period config. Whatever is selected is
#                      transferred to the sandbox under a fixed name
#                      (analysis_config.toml), so runtime_wrapper.sh never needs to
#                      know the source filename — pointing this at a different run
#                      period/collision system's TOML needs no other change here.
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
# Note: the condor Arguments= line generated later in this script (runAsymmetry
# INPUT OUTPUT MODE) stays positional on purpose -- it's an internal,
# script-generated contract consumed by runtime_wrapper.sh, never hand-typed,
# so there's no typo risk to guard against. Only this script's own top-level
# CLI (what a human actually types) needs named arguments.
#
# Prerequisites:
#   - On lxplus, run cmsenv before configuring/building this repo:
#       source /cvmfs/cms.cern.ch/cmsset_default.sh
#       cd <CMSSW_RELEASE>/src && cmsenv && cd /path/to/L2Residuals-2024ppref
#       cmake -S . -B build && cmake --build build
#     or, if the repo is checked out as $CMSSW_BASE/src/Analysis/L2Residuals:
#       scram b -j4
#   - All five cone JEC files must be present in data/jec/ before submitting
#   - The CMSSW src directory used on worker jobs is derived from cmsenv
#     (CMSSW_BASE/src) at submission time, not a TOML setting -- just make
#     sure you're submitting from the cmsenv shell for the release you want.

set -euo pipefail

usage() {
    echo "Usage: $0 -output dir -alltxt -config path [-nosubmit] [-tag value]" >&2
    echo "       $0 -output dir -filelists a.txt b.txt ... -config path [-nosubmit] [-tag value]" >&2
    exit 1
}

if [[ $# -lt 1 ]]; then usage; fi

OUTPUT_DIR=""
USE_ALL=false
NO_SUBMIT=false
EXPLICIT_FILELISTS=()
RUN_TAG=""
CONFIG_PATH=""

while [[ $# -gt 0 ]]; do
    arg="$1"
    if [[ "${arg}" != -?* ]]; then
        echo "ERROR: argument \"${arg}\" is not a -flag (this script only accepts JetMET-style -flag arguments)" >&2
        usage
    fi
    flag="$(printf '%s' "${arg#-}" | tr '[:lower:]' '[:upper:]')"
    shift
    case "${flag}" in
        OUTPUT)
            if [[ $# -eq 0 ]]; then echo "ERROR: -output requires a value" >&2; usage; fi
            OUTPUT_DIR="$1"; shift ;;
        ALLTXT)
            USE_ALL=true ;;
        NOSUBMIT)
            NO_SUBMIT=true ;;
        TAG)
            if [[ $# -eq 0 ]]; then echo "ERROR: -tag requires a value" >&2; usage; fi
            RUN_TAG="$1"; shift ;;
        CONFIG)
            if [[ $# -eq 0 ]]; then echo "ERROR: -config requires a value" >&2; usage; fi
            CONFIG_PATH="$1"; shift ;;
        FILELISTS)
            EXPLICIT_FILELISTS=()
            while [[ $# -gt 0 && "$1" != -?* ]]; do
                EXPLICIT_FILELISTS+=("$1")
                shift
            done
            if [[ ${#EXPLICIT_FILELISTS[@]} -eq 0 ]]; then
                echo "ERROR: -filelists requires at least one value" >&2
                usage
            fi
            ;;
        *)
            echo "ERROR: unrecognized argument \"-${flag}\"" >&2
            usage ;;
    esac
done

if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "ERROR: -output dir is required" >&2
    usage
fi

if [[ "${RUN_TAG}" == */* ]]; then
    echo "ERROR: -tag must not contain '/': ${RUN_TAG}" >&2
    exit 1
fi

if [[ "${USE_ALL}" == true && ${#EXPLICIT_FILELISTS[@]} -gt 0 ]]; then
    echo "ERROR: pass either -alltxt or -filelists, not both" >&2
    usage
fi

if [[ "${USE_ALL}" == false && ${#EXPLICIT_FILELISTS[@]} -eq 0 ]]; then
    echo "ERROR: specify -alltxt or -filelists a.txt b.txt ..." >&2
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
    echo "ERROR: -config path is required -- which TOML gets submitted is a physics-affecting" >&2
    echo "       choice (e.g. main run-period config vs. a closure config), so there is no" >&2
    echo "       implicit default. Pass -config cfg/2024ppRef.toml explicitly if that's what you want." >&2
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

if [[ -z "${CMSSW_BASE:-}" ]]; then
    echo "ERROR: cmsenv is not active. The binary must be built and submitted from a cmsenv shell." >&2
    echo "       source /cvmfs/cms.cern.ch/cmsset_default.sh" >&2
    echo "       cd <CMSSW>/src && cmsenv && cd -" >&2
    echo "       cd /path/to/L2Residuals-2024ppref" >&2
    echo "       cmake -S . -B build && cmake --build build" >&2
    echo "       or build the SCRAM executable with: scram b -j4" >&2
    exit 1
fi

# The CMSSW src directory used for worker jobs is derived from the active
# cmsenv, not a TOML setting -- CMSSW_BASE (set by cmsenv, already required
# above) is the release area, so CMSSW_BASE/src is exactly the src directory
# a worker job needs to `scramv1 runtime -sh` from. This always reflects
# whichever CMSSW release the submitter actually has set up, so there's
# nothing to keep in sync in the TOML and nothing to get stale.
CMSSW_SRC_FROM_CONFIG="${CMSSW_BASE}/src"
if [[ ! -d "${CMSSW_SRC_FROM_CONFIG}" ]]; then
    echo "ERROR: CMSSW_BASE/src does not exist: ${CMSSW_SRC_FROM_CONFIG}" >&2
    exit 1
fi

if [[ ! -f "${BINARY}" ]]; then
    echo "ERROR: no runAsymmetry executable found." >&2
    echo "       For CMake, run: cmake --build build" >&2
    echo "       For SCRAM, put the repo at $CMSSW_BASE/src/Analysis/L2Residuals and run: scram b -j4" >&2
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
    FILELISTS=()
    for f in "${EXPLICIT_FILELISTS[@]}"; do
        if [[ -d "${f}" ]]; then
            DIR_TXT=("${f}"/*.txt)
            if [[ ! -f "${DIR_TXT[0]}" ]]; then
                echo "ERROR: no .txt filelists found in directory: ${f}" >&2
                exit 1
            fi
            FILELISTS+=("${DIR_TXT[@]}")
        elif [[ -f "${f}" ]]; then
            FILELISTS+=("${f}")
        else
            echo "ERROR: filelist not found (not a file or directory): ${f}" >&2
            exit 1
        fi
    done
fi

# Resolve to absolute paths now, before cd-ing into WORKDIR below --
# filelists are read directly from their original location (never copied
# into the submission sandbox; see the data/ staging comment further down
# for why data/txt/ specifically isn't part of that copy).
for i in "${!FILELISTS[@]}"; do
    FILELISTS[$i]="$(cd "$(dirname "${FILELISTS[$i]}")" && pwd)/$(basename "${FILELISTS[$i]}")"
done

TODAY=$(date +"%Y-%m-%d_%H-%M-%S")
SUBMISSIONS_DIR="${CONDOR_DIR}/submissions"
mkdir -p "${SUBMISSIONS_DIR}"
WORKDIR="${SUBMISSIONS_DIR}/${TODAY}"
mkdir -p "${WORKDIR}"
# asymmetry_<TAG> keeps separate passes (e.g. -tag abs_eta vs -tag clos_dir_eta)
# from landing in the same output tree; plain "asymmetry" when no -tag is given.
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
    # Only jec/veto/json are read by the worker (via [paths]/[jec] in
    # analysis_config.toml, resolved relative to this sandbox's data/) --
    # data/txt/ is read once here on the submit host to build the per-job
    # Arguments= lines below, never by the worker itself, and data/root/
    # and data/sample/ are local-only derived/sample artifacts a worker
    # never touches at all. Copying the whole data/ directory here used to
    # stage several GB of dead weight into every submission sandbox.
    mkdir -p data
    for sub in jec veto json; do
        if [[ -d "${DATA_DIR}/${sub}" ]]; then
            cp -r "${DATA_DIR}/${sub}" "data/${sub}"
        fi
    done
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

        FILELIST_FILE="${FILELIST_PATH}"
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

Transfer_Input_Files    = $(pwd)/runtime_wrapper.sh,$(pwd)/runAsymmetry${LIBRARY:+,$(pwd)/libl2residuals.so},$(pwd)/data,$(pwd)/analysis_config.toml

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

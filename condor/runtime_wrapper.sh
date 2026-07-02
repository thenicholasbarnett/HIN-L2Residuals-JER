#!/bin/bash -l
# Executed on each Condor worker node.
# Sets up the CMSSW ROOT environment, then runs runAsymmetry.
#
# Called by Condor as:
#   runtime_wrapper.sh runAsymmetry <input.root> <output.root> <mode>

set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 EXECUTABLE INPUT OUTPUT MODE" >&2
    exit 1
fi

EXECUTABLE="$1"
INPUT="$2"
OUTPUT="$3"
MODE="$4"
START_DIR="$(pwd)"

CMSSW_SRC="@CMSSW_SRC@"

if [[ -z "${CMSSW_SRC}" || "${CMSSW_SRC}" == "@CMSSW_SRC@" ]]; then
    echo "ERROR: CMSSW_SRC was not stamped into runtime_wrapper.sh by make_condor.sh" >&2
    exit 1
fi

echo "CMSSW environment: $(basename "$(dirname "${CMSSW_SRC}")")"
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${CMSSW_SRC}"
eval "$(scramv1 runtime -sh)"
cd "${START_DIR}"
export L2RESIDUALS_CONFIG="${START_DIR}/analysis_config.toml"

# Binary is compiled against CMSSW ROOT (which has libCling).
# Prepend . for libl2residuals.so transferred to the sandbox.
export LD_LIBRARY_PATH="./${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

echo "Input:  ${INPUT}"
echo "Output: ${OUTPUT}"
echo "Mode:   ${MODE}"

chmod +x "${EXECUTABLE}"
# The compiled binary's CLI is entirely KEY=value tokens, no positional args.
./"${EXECUTABLE}" "INPUT=${INPUT}" "OUTPUT=${OUTPUT}" "MODE=${MODE}" "CONFIG=${START_DIR}/analysis_config.toml"

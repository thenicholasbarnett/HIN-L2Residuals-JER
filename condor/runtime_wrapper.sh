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

# Set this to your CMSSW src directory on AFS, e.g.:
# /afs/cern.ch/user/n/nbarnett/public/condor/workArea/CMSSW_13_2_4/src
CMSSW_SRC=""

if [[ -z "${CMSSW_SRC}" ]]; then
    echo "ERROR: CMSSW_SRC is not set in condor/runtime_wrapper.sh" >&2
    exit 1
fi

echo "CMSSW environment: $(basename "$(dirname "${CMSSW_SRC}")")"
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${CMSSW_SRC}"
eval "$(scramv1 runtime -sh)"
cd "${START_DIR}"

# The shared library libl2residuals.so was transferred to the sandbox directory.
export LD_LIBRARY_PATH=.:${LD_LIBRARY_PATH:-}

echo "Input:  ${INPUT}"
echo "Output: ${OUTPUT}"
echo "Mode:   ${MODE}"

chmod +x "${EXECUTABLE}"
./"${EXECUTABLE}" "${INPUT}" "${OUTPUT}" "${MODE}"

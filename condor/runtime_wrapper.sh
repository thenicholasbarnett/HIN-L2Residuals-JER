#!/bin/bash -l
# Executed on each Condor worker node
# Sets up the CMSSW ROOT environment, then runs runAsymmetry.

set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 EXECUTABLE INPUT OUTPUT MODE CMSSW_SRC" >&2
  exit 1
fi

EXECUTABLE="$1"
INPUT="$2"
OUTPUT="$3"
MODE="$4"
CMSSW_SRC="$5"
START_DIR="$(pwd)"

if [[ -z "${CMSSW_SRC}" ]]; then
  echo "ERROR: CMSSW_SRC argument was empty" >&2
  exit 1
fi

echo "CMSSW environment: $(basename "$(dirname "${CMSSW_SRC}")")"
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${CMSSW_SRC}"
eval "$(scramv1 runtime -sh)"
cd "${START_DIR}"
export L2RESIDUALS_CONFIG="${START_DIR}/analysis_config.toml"

# Binary is compiled against CMSSW ROOT
export LD_LIBRARY_PATH="./${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

echo "Input:  ${INPUT}"
echo "Output: ${OUTPUT}"
echo "Mode:   ${MODE}"

chmod +x "${EXECUTABLE}"
# "-key value" CommandLine parser
./"${EXECUTABLE}" -input "${INPUT}" -output "${OUTPUT}" -mode "${MODE}" -config "${START_DIR}/analysis_config.toml"

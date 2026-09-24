#!/bin/bash -l
# Executed on each Condor worker node
# Sets up the CMSSW ROOT environment, then runs runAsymmetry and/or
# runDijetProfiles on the same input (an output of "none" skips that binary).

set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: $0 INPUT MODE CMSSW_SRC CLOSURE ASYM_OUTPUT PROFILES_OUTPUT" >&2
  exit 1
fi

INPUT="$1"
MODE="$2"
CMSSW_SRC="$3"
CLOSURE="$4"
ASYM_OUTPUT="$5"
PROFILES_OUTPUT="$6"
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
echo "Asymmetry output: ${ASYM_OUTPUT}"
echo "Profiles output:  ${PROFILES_OUTPUT}"
echo "Mode:   ${MODE}"
echo "Closure: ${CLOSURE}"

# "-key value" CommandLine parser
CLOSURE_ARGS=()
if [[ "${CLOSURE}" == "true" ]]; then
  CLOSURE_ARGS=(-calibration jer -closure true)
fi
if [[ "${ASYM_OUTPUT}" != none ]]; then
  chmod +x runAsymmetry
  ./runAsymmetry -input "${INPUT}" -output "${ASYM_OUTPUT}" -mode "${MODE}" -config "${START_DIR}/analysis_config.toml" "${CLOSURE_ARGS[@]}"
fi
if [[ "${PROFILES_OUTPUT}" != none ]]; then
  chmod +x runDijetProfiles
  ./runDijetProfiles -input "${INPUT}" -output "${PROFILES_OUTPUT}" -mode "${MODE}" -config "${START_DIR}/analysis_config.toml" "${CLOSURE_ARGS[@]}"
fi

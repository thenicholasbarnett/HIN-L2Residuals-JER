#!/usr/bin/env bash
# Verifies the compiled executable and interpreted macro are both usable.
# Run from the repo root after: cmake --build build

set -euo pipefail

PASS=0
FAIL=0

check() {
    # Using PASS=$((PASS+1)) instead of ((PASS++)) to avoid set -e triggering
    # when the pre-increment value is 0 (bash arithmetic returns exit 1 for 0).
    if eval "$2" &>/dev/null; then
        echo "  PASS  $1"
        PASS=$((PASS + 1))
    else
        echo "  FAIL  $1"
        FAIL=$((FAIL + 1))
    fi
}

# Shared library has .so on Linux and .dylib on macOS
case "$(uname)" in
    Darwin) LIB_EXT="dylib" ;;
    *)      LIB_EXT="so"    ;;
esac
LIB="lib/libl2residuals.${LIB_EXT}"

echo ""
echo "=== TestBuild ==="

echo ""
echo "[compiled]"
check "executable exists"            "test -f bin/runAsymmetry"
check "executable is runnable"       "test -x bin/runAsymmetry"
check "prints usage with no args"    "{ bin/runAsymmetry 2>&1 || true; } | grep -q 'Usage:'"
check "exits non-zero with no args"  "! bin/runAsymmetry"
check "runTextFile exists"           "test -f bin/runTextFile"
check "runResiduals exists"          "test -f bin/runResiduals"
check "plotResiduals exists"         "test -f bin/plotResiduals"

echo ""
echo "[interpreted]"
check "library exists"               "test -f ${LIB}"
check "ROOT can load library"        "root -l -b -q -e 'gSystem->Load(\"${LIB}\"); gApplication->Terminate(0)'"
check "ROOT can load macro"          "root -l -b -q -e '.include include' -e '.include cfg' -e 'gSystem->Load(\"${LIB}\")' -e '.L macros/runAsymmetry.C'"

echo ""
echo "=== ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]

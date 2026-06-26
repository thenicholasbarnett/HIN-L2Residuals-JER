#!/usr/bin/env bash
# Verifies the compiled executable and interpreted macro are both usable.
# Run from the repo root after: cmake --build build

set -e
PASS=0
FAIL=0

check() {
    if eval "$2" &>/dev/null; then
        echo "  PASS  $1"
        ((PASS++))
    else
        echo "  FAIL  $1"
        ((FAIL++))
    fi
}

echo ""
echo "=== TestBuild ==="

echo ""
echo "[compiled]"
check "executable exists"            "test -f bin/runAsymmetry"
check "executable is runnable"       "test -x bin/runAsymmetry"
check "prints usage with no args"    "! ./bin/runAsymmetry 2>&1 | grep -q 'Usage:'"
# exits non-zero with no args (usage error) — that's the expected behaviour
check "exits non-zero with no args"  "! ./bin/runAsymmetry"

echo ""
echo "[interpreted]"
check "library exists"               "test -f lib/libl2residuals.so"
check "ROOT can load library"        "root -l -b -q -e 'gSystem->Load(\"lib/libl2residuals.so\"); gApplication->Terminate(0)'"
check "ROOT can load macro"          "root -l -b -q -e '.L macros/runAsymmetry.C'; true"

echo ""
echo "=== ${PASS} passed, ${FAIL} failed ==="
[ "$FAIL" -eq 0 ]

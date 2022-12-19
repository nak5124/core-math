#!/bin/bash

set -e

if [ -z "$LAST_COMMIT" ]; then
    LAST_COMMIT="HEAD~"
    if [ -n "$CI_COMMIT_BEFORE_SHA" ] && [ "$CI_COMMIT_BEFORE_SHA" != "0000000000000000000000000000000000000000" ]; then
        LAST_COMMIT="$CI_COMMIT_BEFORE_SHA"
    fi
fi

FUNCTIONS_EXHAUSTIVE=(acosf acoshf acospif asinf asinhf asinpif atanf atanhf atanpif cbrtf cosf coshf cospif erff erfcf expf exp10f exp10m1f exp2f exp2m1f expm1f logf log10f log10p1f log1pf log2f log2p1f rsqrtf sinf sinhf sinpif tanf tanhf tanpif)
FUNCTIONS_WORST=(acos asin atan2f cbrt exp exp2 hypotf log rsqrt)
FUNCTIONS_SPECIAL=(hypotf)

echo "Reference commit is $LAST_COMMIT"

check () {
    KIND="$1"
    if ! { echo "$FORCE_FUNCTIONS" | tr ' ' '\n' | grep --quiet '^'"$FUNCTION"'$'; } && git diff --quiet "$LAST_COMMIT".. -- src/*/*/$FUNCTION.c; then
        echo "Skipped $FUNCTION"
    else
        echo "Checking $FUNCTION..."
        ./check.sh "$KIND" "$FUNCTION"
    fi
}

for FUNCTION in "${FUNCTIONS_EXHAUSTIVE[@]}"; do
    check --exhaustive
done

for FUNCTION in "${FUNCTIONS_WORST[@]}"; do
    check --worst
done

for FUNCTION in "${FUNCTIONS_SPECIAL[@]}"; do
    check --special
done

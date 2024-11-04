#!/bin/bash
# Examples of usage:
# DRY=--dry ./ci.sh to only try compilation (of last modified functions)
# FORCE=true DRY=--dry ./ci.sh to only try compilation (of all functions)
# FORCE_FUNCTIONS="xxx yyy" ./ci.sh to force checking xxx and yyy
# CC=clang CFLAGS=-Werror ./ci.sh

set -e # We want the above command to potentially fail, only set -e now.

if [ -z "$LAST_COMMIT" ]; then
    LAST_COMMIT="HEAD~"
    if [ -n "$CI_COMMIT_BEFORE_SHA" ] && [ "$CI_COMMIT_BEFORE_SHA" != "0000000000000000000000000000000000000000" ]; then
        LAST_COMMIT="$CI_COMMIT_BEFORE_SHA"
    fi
fi

# use the same order as on https://core-math.gitlabpages.inria.fr/
FUNCTIONS_EXHAUSTIVE=(acosf acoshf acospif asinf asinhf asinpif atanf atanhf atanpif cbrtf cosf coshf cospif erff erfcf expf exp10f exp10m1f exp2f exp2m1f expm1f lgammaf logf log10f log10p1f log1pf log2f log2p1f rsqrtf sincosf sinf sinhf sinpif tanf tanhf tanpif tgammaf)
FUNCTIONS_WORST=(acos acosh acospi asin asinh asinpi atan atan2 atan2f atan2pi atan2pif atanh atanpi cbrt cbrtl cos cosh cospi erf erfc exp expl exp10 exp10m1 exp2 exp2l exp2m1 hypot hypotf log log10 log10p1 log1p log2 log2l log2p1 pow powf powl rsqrt rsqrtl sin sincos sinh sinpi tan tanh tanpi)
FUNCTIONS_SPECIAL=(acos acosh acospi asin asinh asinpi atan atan2 atan2f atan2pi atan2pif atanh atanpi cbrt cos cosh cospi erf erfc exp exp10 exp10m1 exp10m1f exp2 exp2m1 exp2m1f expm1 hypot hypotf log log10 log10p1 log1p log2 log2p1 pow powf powl rsqrt sin sinh sinpi tan tanh tanpi)

echo "Reference commit is $LAST_COMMIT"

check () {
    KIND="$1"
    if [ "$FORCE" != "" ]; then
        doit=1
    elif ! { echo "$FORCE_FUNCTIONS" | tr ' ' '\n' | grep --quiet '^'"$FUNCTION"'$'; } && git diff --quiet "$LAST_COMMIT".. -- src/*/*/$FUNCTION.c; then
        doit=0
    else
        doit=1
    fi
		if [ "$doit" == "1" ] && [ "$SKIP128" == "1" ] && $CC -E src/*/*/$FUNCTION.c | grep -q  __int128; then
        echo "__int128 support is needed for" $FUNCTION "but is not available"
        doit=0
    fi
    if [ "$doit" == "0" ]; then
        echo "Skip $FUNCTION"
    else
        echo "Checking $FUNCTION..."				
	# we want to detect compiler warnings
        EXTRA_CFLAGS=-Werror ./check.sh $DRY "$KIND" "$FUNCTION"
    fi
}

if [ -z "$CC" ]; then
	CC="cc"
fi

if $CC -E $CFLAGS ci/int128test.c -o /dev/null &> /dev/null; then
    echo "Compiler supports __int128"
else
    echo "Compiler lacks __int128 support"
    SKIP128=1
fi

for FUNCTION in "${FUNCTIONS_EXHAUSTIVE[@]}"; do
    check --exhaustive
done

for FUNCTION in "${FUNCTIONS_WORST[@]}"; do
    check --worst
done

for FUNCTION in "${FUNCTIONS_SPECIAL[@]}"; do
    check --special
done

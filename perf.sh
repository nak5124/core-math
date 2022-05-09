#!/bin/bash
# Usage: ./perf.sh acos

RANDOMS_FILE="$(mktemp /tmp/core-math.XXXXXX)"
trap "rm -f $RANDOMS_FILE" 0

f=$1
u="$(echo src/binary*/*/$f.c)"

if [ -z "$CORE_MATH_PERF_MODE" ]; then
    echo 'Please set CORE_MATH_PERF_MODE environment variable (perf or rdtsc)'
    exit 2
fi

if [ -f "$u" ]; then
    dir="${u%/*}"
else
    echo "Unknown function: $f"
    exit 2
fi

# for clang we might want to add -ffp-contract=on to enable FMA
# -ffinite-math-only is needed to inline fmaxf and fminf
if [ "$CFLAGS" == "" ]; then
   export CFLAGS="-O3 -march=native -ffinite-math-only"
fi

if [ -n "$LIBM" ]; then
    BACKUP_LIBM="$LIBM"
    unset LIBM
fi

if [ -z "$CORE_MATH_QUIET" ]; then
    make -s -C src/generic/support clean
    make -s -C src/generic/support all
    $CORE_MATH_LAUNCHER src/generic/support/glibc_version >&2
fi

cd $dir
make -s clean
make -s perf
./perf --file "$RANDOMS_FILE" --reference --count 1000000

if [ "$CORE_MATH_PERF_MODE" = perf ]; then
    perf stat -e cpu-cycles --no-big-num -x';' $CORE_MATH_LAUNCHER ./perf --file "$RANDOMS_FILE" --count 1000000 --repeat 1000 $PERF_ARGS 2>&1 | { IFS=';' read a b; printf 'scale=3\n%s/1000000000\n' $a | bc; }
    perf stat -e cpu-cycles --no-big-num -x';' $CORE_MATH_LAUNCHER ./perf --file "$RANDOMS_FILE" --count 1000000 --repeat 1000 --libc $PERF_ARGS 2>&1 | { IFS=';' read a b; printf 'scale=3\n%s/1000000000\n' $a | bc; }
elif [ "$CORE_MATH_PERF_MODE" = rdtsc ]; then
    $CORE_MATH_LAUNCHER ./perf --file "$RANDOMS_FILE" --count 1000000 --repeat 1000 $PERF_ARGS --rdtsc
    $CORE_MATH_LAUNCHER ./perf --file "$RANDOMS_FILE" --count 1000000 --repeat 1000 --libc $PERF_ARGS --rdtsc
fi

if [ -n "$BACKUP_LIBM" ]; then
    export LIBM="$BACKUP_LIBM"
    make -s clean
    make -s perf
    if [ "$CORE_MATH_PERF_MODE" = perf ]; then
        perf stat -e cpu-cycles --no-big-num -x';' $CORE_MATH_LAUNCHER ./perf --file "$RANDOMS_FILE" --count 1000000 --repeat 1000 --libc $PERF_ARGS 2>&1 | { IFS=';' read a b; printf 'scale=3\n%s/1000000000\n' $a | bc; }
    elif [ "$CORE_MATH_PERF_MODE" = rdtsc ]; then
        $CORE_MATH_LAUNCHER ./perf --file "$RANDOMS_FILE" --count 1000000 --repeat 1000 --libc $PERF_ARGS --rdtsc
    fi
fi

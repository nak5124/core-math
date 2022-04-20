#!/bin/bash
# Usage: ./perf.sh acos

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

cd $dir
make -s clean
make -s perf
./perf --file /tmp/randoms.dat --reference --count 1000000

if [ "$CORE_MATH_PERF_MODE" = perf ]; then
    perf stat -e cpu-cycles --no-big-num -x';' ./perf --file /tmp/randoms.dat --count 1000000 --repeat 1000 $PERF_ARGS 2>&1 | { IFS=';' read a b; printf 'scale=3\n%s/1000000000\n' $a | bc; }
    perf stat -e cpu-cycles --no-big-num -x';' ./perf --file /tmp/randoms.dat --count 1000000 --repeat 1000 --libc $PERF_ARGS 2>&1 | { IFS=';' read a b; printf 'scale=3\n%s/1000000000\n' $a | bc; }
elif [ "$CORE_MATH_PERF_MODE" = rdtsc ]; then
    ./perf --file /tmp/randoms.dat --count 1000000 --repeat 1000 $PERF_ARGS --rdtsc
    ./perf --file /tmp/randoms.dat --count 1000000 --repeat 1000 --libc $PERF_ARGS --rdtsc
fi

if [ -n "$BACKUP_LIBM" ]; then
    export LIBM="$BACKUP_LIBM"
    make -s clean
    make -s perf
    if [ "$CORE_MATH_PERF_MODE" = perf ]; then
        perf stat -e cpu-cycles --no-big-num -x';' ./perf --file /tmp/randoms.dat --count 1000000 --repeat 1000 --libc $PERF_ARGS 2>&1 | { IFS=';' read a b; printf 'scale=3\n%s/1000000000\n' $a | bc; }
    elif [ "$CORE_MATH_PERF_MODE" = rdtsc ]; then
        ./perf --file /tmp/randoms.dat --count 1000000 --repeat 1000 --libc $PERF_ARGS --rdtsc
    fi
fi

#!/bin/bash
# Usage: ./perf.sh acos

f=$1
u="$(echo src/binary*/*/$f.c)"

if [ -f "$u" ]; then
    dir="${u%/*}"
else
    echo "Unknown function: $f"
    exit 2
fi

export CFLAGS="-O3 -march=native"

if [ -n "$LIBM" ]; then
    BACKUP_LIBM="$LIBM"
    unset LIBM
fi

cd $dir
make -s clean
make -s perf
./perf --file /tmp/randoms.dat --reference --count 1000000
perf stat -e cpu-cycles --no-big-num -x';' ./perf --file /tmp/randoms.dat --count 1000000 --repeat 1000 $PERF_ARGS 2>&1 | { IFS=';' read a b; printf 'scale=3\n%s/1000000000\n' $a | bc; }
perf stat -e cpu-cycles --no-big-num -x';' ./perf --file /tmp/randoms.dat --count 1000000 --repeat 1000 --libc $PERF_ARGS 2>&1 | { IFS=';' read a b; printf 'scale=3\n%s/1000000000\n' $a | bc; }

if [ -n "$BACKUP_LIBM" ]; then
    export LIBM="$BACKUP_LIBM"
    make -s clean
    make -s perf
    perf stat -e cpu-cycles --no-big-num -x';' ./perf --file /tmp/randoms.dat --count 1000000 --repeat 1000 --libc $PERF_ARGS 2>&1 | { IFS=';' read a b; printf 'scale=3\n%s/1000000000\n' $a | bc; }
fi

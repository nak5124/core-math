#!/bin/bash
# Usage: ./throughput.sh acos
f=$1
cd src/binary32/$f
make -s clean
make -s CFLAGS="-O3 -march=native"
./perf --file /tmp/randoms.dat --reference --count 1000000
perf stat -e cpu-cycles --no-big-num -x';' ./perf --file /tmp/randoms.dat --count 1000000 --repeat 1000 2>&1 | { IFS=';' read a b; printf 'scale=3\n%s/1000000000\n' $a | bc; }

#!/bin/bash
# Usage: ./perf.sh acos

S=20 # trial
N=100000 # count
M=500 # repeat

read -r -d '' prog_end <<EOF
END {
  s = 0;
  nout = int(5*i/100);
  for(k = 1; k < i-nout; k++){
    d = a[k] - a[0];
    s += d*d;
  }
  rms = sqrt(s/(i-nout-1));
  nmd = a[int(i/2)] - a[0];
  printf "Ntrial = %d ; Min = %.3f + %.3f clc/call; Median-Min = %.3f clc/call; Max = %.3f clc/call;\n", i, a[0], rms, nmd, a[i-1];
}
EOF

prog_bar () {
    local T=$1 i=$2
    local p=$(( i * 100 / T ))
    local j=$(( T - i ))
    printf "[%s%s] %3d %%\r" "${str_hsh:0:${i}}" "${str_dot:0:${j}}" $p
    if (( i == T )); then
       echo
    fi
}

collect_perf_stat () {
    echo -n "" > $LOG_FILE
    if [ -z "$CORE_MATH_QUIET" ]; then
	prog_bar $S 0
    fi
    local i=1
    while [ $i -le $S ]; do
	perf stat -e cycles -x " " ./perf $PERF_ARGS &>> $LOG_FILE
	if [ -z "$CORE_MATH_QUIET" ]; then
	    prog_bar $S $i
	fi
	i=$(( i + 1 ))
    done
}
process_perf_stat () {
    sort -g -k 1 $LOG_FILE | awk "/cycles:u/{a[i++]=\$1/(${N}*${M});} ${prog_end}"
}
proc_perf () {
    collect_perf_stat
    process_perf_stat
}

collect_rdtsc_stat () {
    echo -n "" > $LOG_FILE
    local i=1
    while [ $i -le $S ]; do
	./perf $PERF_ARGS &>> $LOG_FILE
	i=$(( i + 1 ))
    done
}
process_rdtsc_stat () {
    sort -g $LOG_FILE | awk "{a[i++]=\$1;} ${prog_end}"
}
proc_rdtsc () {
    collect_rdtsc_stat
    process_rdtsc_stat
}

has_symbol () {
    [ "$(nm "$LIBM" | while read a b c; do if [ "$c" = "$f" ]; then echo OK; return; fi; done | wc -l)" -ge 1 ]
}

RANDOMS_FILE="$(mktemp /tmp/core-math.XXXXXX)"
LOG_FILE="$(mktemp /tmp/core-math.XXXXXX)"
trap "rm -f $RANDOMS_FILE $LOG_FILE" 0

f=$1
u="$(echo src/binary*/*/$f.c)"

if [ -z "$CORE_MATH_PERF_MODE" ]; then
    echo 'CORE_MATH_PERF_MODE (perf or rdtsc) environment variable is not set. The default is perf.'
    CORE_MATH_PERF_MODE=perf
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
    str_hsh=""
    str_dot=""
    i=1
    while [ $i -le $S ]; do
	str_hsh="${str_hsh}#"
	str_dot="${str_dot}."
	i=$(( i + 1 ))
    done
fi

cd $dir
make -s clean
make -s perf

# prepare random arguments for performance test
./perf --file ${RANDOMS_FILE} --count ${N} --reference

PERF_ARGS="${PERF_ARGS} --file ${RANDOMS_FILE} --count ${N} --repeat ${M}"

if [ "$CORE_MATH_PERF_MODE" = perf ]; then
    proc_perf

    PERF_ARGS="${PERF_ARGS} --libc"
    proc_perf

elif [ "$CORE_MATH_PERF_MODE" = rdtsc ]; then
    PERF_ARGS="${PERF_ARGS} --rdtsc"
    proc_rdtsc

    PERF_ARGS="${PERF_ARGS} --libc"
    proc_rdtsc
fi

if [ -n "$BACKUP_LIBM" ]; then
    export LIBM="$BACKUP_LIBM"
    if has_symbol; then
	PERF_ARGS="${PERF_ARGS} --libc"
	make -s clean
	make -s perf
	if [ "$CORE_MATH_PERF_MODE" = perf ]; then
	    proc_perf
	    
	elif [ "$CORE_MATH_PERF_MODE" = rdtsc ]; then
	    PERF_ARGS="${PERF_ARGS} --rdtsc"
	    proc_rdtsc
	fi
    elif [ -z "$CORE_MATH_QUIET" ]; then
        echo "$f is not present in $LIBM; skipping" >&2
    fi
fi

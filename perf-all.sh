#!/bin/bash

if [ -z "$CORE_MATH_QUIET" ]; then
    make -s -C src/generic/support clean
    make -s -C src/generic/support all
    $CORE_MATH_LAUNCHER src/generic/support/glibc_version >&2
fi

for u in src/binary*/*/Makefile; do
    f="$(sed -n 's/FUNCTION_UNDER_TEST := //p' $u)"
    echo -n "$f "
    CORE_MATH_QUIET=1 ./perf.sh $f | xargs echo
done

#!/bin/bash

for u in src/binary*/*/Makefile; do
    f="$(sed -n 's/FUNCTION_UNDER_TEST := //p' $u)"
    echo -n "$f "
    ./perf.sh $f | xargs echo
done

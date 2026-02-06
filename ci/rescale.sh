#!/bin/bash

for f in acos acosh acospi asin asinh asinpi atan atanh atanpi cbrt cos cosh cospi erf erfc exp exp10 exp10m1 exp2 exp2m1 expm1 lgamma log log10 log10p1 log1p log2 log2p1 rsqrt sin sinh sinpi support tan tanh tanpi; do
   e=`shuf -i 0-2098 -n 1`
   e=`echo $e - 1074 | bc -q`
   echo Testing $f with rescale $e
   ./check.sh --worst --rescale $e $f
done

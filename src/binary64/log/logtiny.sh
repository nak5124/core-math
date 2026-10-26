#!/bin/bash
f=`mktemp`
for r in ne zr up dn; do
   echo $r
   sed "s/RND/$r/g" logtiny.g > $f
   gappa $f
done
rm $f

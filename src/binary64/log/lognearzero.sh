#!/bin/bash
f=`mktemp`
for r in ne zr up dn; do
   sed "s/RND/$r/g" lognearzero.g | sed "s/XMIN/0x1.0000000000001p+0/g" | sed "s/XMAX/0x1.008p+0/g" > $f
   gappa $f
   sed "s/RND/$r/g" lognearzero.g | sed "s/XMIN/0x1.ffp-1/g" | sed "s/XMAX/0x1.fffffffffffffp-1/g" > $f
   gappa $f
done
rm $f

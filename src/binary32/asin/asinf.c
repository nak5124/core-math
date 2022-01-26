/* Correctly-rounded arc-sine function for binary32 value.

Copyright (c) 2022 Alexei Sibidanov.

This file is part of the CORE-MATH project
(https://core-math.gitlabpages.inria.fr/).

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Tested on x86_64-linux with and without FMA (-march=native). */

/* Notes:
 * this function is supposed to deal properly with special arguments as well
   as raise proper exceptions, but this is not fully tested.
 * the 'main' code tests exhaustively all binary32 values, with -rndn by
   default (rounding to nearest), -rndz (towards zero), -rndu (upwards)
   or -rndd (downwards);
 * tested on x86-64-linux with gcc 10.2.1, with and without -march=native.
*/

#include <stdint.h>
#include <errno.h>
#include <fenv.h>

typedef union {float f; unsigned u;} b32u32_u;

float cr_asinf(float x){
  const double pi = 0x1.921fb54442d18p+1;
  b32u32_u t = {.f = x};
  int e = (t.u>>23)&0xff;
  double r;
  if(__builtin_expect(e>=127, 0)){
    if(__builtin_fabsf(x)==1.0f) return __builtin_copysign(pi/2, (double)x);
    if(e==0xff && (t.u<<9)) return x; // nan
    errno = EDOM;
    feraiseexcept(FE_INVALID);
    return __builtin_nanf("1");
  }
  if (e<126){
    if (__builtin_expect(e<127-12, 0))
      return __builtin_fmaf(x, 0x1p-25, x);
    static const double c[] =
      {0x1.555555555529cp-3, 0x1.333333337e0ddp-4, 0x1.6db6db3b4465ep-5, 0x1.f1c72e13ac306p-6,
       0x1.6e89cebe06bc4p-6, 0x1.1c6dcf5289094p-6, 0x1.c6dbbcc7c6315p-7, 0x1.8f8dc2615e996p-7,
       0x1.a5833b7bf15e8p-8, 0x1.43f44ace1665cp-6, -0x1.0fb17df881c73p-6, 0x1.07520c026b2d6p-5};
    double z = x, z2 = z*z, z4 = z2*z2, z8 = z4*z4;
    double c0 = c[0] + z2*c[1];
    double c2 = c[2] + z2*c[3];
    double c4 = c[4] + z2*c[5];
    double c6 = c[6] + z2*c[7];
    double c8 = c[8] + z2*c[9];
    double c10 = c[10] + z2*c[11];
    c0 += c2*z4;
    c4 += c6*z4;
    c8 += c10*z4;
    c0 += z8*(c4 + z8*c8);
    r = z + (z*z2)*c0;
  } else {
    static const struct {union{float arg; unsigned uarg;}; float rh, rl;} st[] = {
      {{0x1.107434p-1f}, 0x1.1f4b64p-1f, 0x1.fffffep-26f},
      {{0x1.55688ap-1f}, 0x1.75b8a2p-1f, 0x1.c5ca9ap-53f}
    };
    unsigned at = t.u&(~0u>>1), sgn = t.u>>31;
    for(int i=0;i<2;i++) {
      if(__builtin_expect(st[i].uarg == at, 0)){
	if(sgn)
	  return -st[i].rh - st[i].rl;
	else
	  return  st[i].rh + st[i].rl;
      }
    }
    float ax = __builtin_fabsf(x);
    double z = 1.0 - ax, s = __builtin_sqrt(z);
    static const double c[] =
      {0x1.6a09e667f3bcbp+0, 0x1.e2b7dddff2db9p-4, 0x1.b27247ab42dbcp-6, 0x1.02995cc4e0744p-7,
       0x1.5ffb0276ec8eap-9, 0x1.033885a928decp-10, 0x1.911f2be23f8c7p-12, 0x1.4c3c55d2437fdp-13,
       0x1.af477e1d7b461p-15, 0x1.abd6bdff67dcbp-15, -0x1.1717e86d0fa28p-16, 0x1.6ff526de46023p-16};
    double z2 = z*z, z4 = z2*z2;
    double c0 = c[0] + z*c[1];
    double c2 = c[2] + z*c[3];
    double c4 = c[4] + z*c[5];
    double c6 = c[6] + z*c[7];
    double c8 = c[8] + z*c[9];
    double c10 = c[10] + z*c[11];
    c0 += c2*z2;
    c4 += c6*z2;
    c8 += z2*c10;
    c0 += z4*(c4 + z4*c8);
    r = pi/2 - s*c0;
    r = __builtin_copysign(r,(double)x);
  }
  return r;
}

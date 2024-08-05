/* Correctly-rounded Euclidean distance function (hypot) for binary32 values.

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
*/

#include <stdint.h>
#include <errno.h>

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union {float f; uint32_t u;} b32u32_u;
typedef union {double f; uint64_t u;} b64u64_u;

float cr_hypotf(float x, float y){
  float ax = __builtin_fabsf(x), ay = __builtin_fabsf(y);
  b32u32_u tx = {.f = ax}, ty = {.f = ay};
  if(__builtin_expect(tx.u >= (0xffu<<23) || ty.u >= (0xffu<<23), 0)){
    if(tx.u == (0xffu<<23)) return ax;
    if(ty.u == (0xffu<<23)) return ay;
    return ax + ay;
  }
  float at = __builtin_fmaxf(ax, ay);
  ay = __builtin_fminf(ax, ay);
  double xd = at, yd = ay, x2 = xd*xd, y2 = yd*yd, r2 = x2 + y2;
  if(__builtin_expect(yd < xd*0x1.fffffep-13, 0)) return __builtin_fmaf(0x1p-13f, ay, at);
  double r = __builtin_sqrt(r2);
  b64u64_u t = {.f = r};
  float c = r;
  if(t.u>(uint64_t)0x47efffffe0000000ul){
    if(c>0x1.fffffep127f) errno = ERANGE;
    return c;
  }
  if(__builtin_expect(((t.u + 1)&0xfffffff) > 2, 1)) return c;
  double cd = c;
  if((cd*cd - x2) - y2 == 0.0) return c;
  double ir2 = 0.5/r2, dr2 = (x2 - r2) + y2;
  double rs = r*ir2, dz = dr2 - __builtin_fma(r, r, -r2), dr = rs*dz;
  double rh = r + dr, rl = dr + (r - rh);
  t.f = rh;
  if(__builtin_expect((t.u&0xfffffff) == 0, 1)){
    if(rl>0.0) t.u += 1;
    if(rl<0.0) t.u -= 1;
  }
  return t.f;
}

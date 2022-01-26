/* Correctly-rounded arc-tangent of binary32 value.

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

Tested on x86_64-linux with and without FMA (-march=native).
*/

#include <stdint.h>

typedef union {float f; unsigned u;} b32u32_u;
typedef union {double f; unsigned long u;} b64u64_u;
typedef unsigned long u64;

float cr_atanf(float x){
  const double pi2 = 0x1.921fb54442d18p+0;
  b32u32_u t = {.f = x};
  int e = (t.u>>23)&0xff, gt = e>=127, sgn = t.u>>31;
  if(__builtin_expect(e==0xff, 0)) {
    if(t.u<<9) return x; // nan
    return __builtin_copysign(pi2,(double)x); // inf
  }
  if (__builtin_expect(e<127-13, 0)){
    if (__builtin_expect(e<127-25, 0))
      return __builtin_fmaf(-x, __builtin_fabsf(x), x);
    return __builtin_fmaf(-0x1.5555555555555p-2f*x, x*x, x);
  }
  unsigned ax = t.u&(~0u>>1);
  double z = x;
  if(gt) z = 1/z;
  double z2 = z*z, z4 = z2*z2, z8 = z4*z4;
  static const double cn[] =
    {0x1p+0, 0x1.40e0698f94c35p+1, 0x1.248c5da347f0dp+1, 0x1.d873386572976p-1, 0x1.46fa40b20f1dp-3,
     0x1.33f5e041eed0fp-7, 0x1.546bbf28667c5p-14};
  static const double cd[] =
    {0x1p+0, 0x1.6b8b143a3f6dap+1, 0x1.8421201d18ed5p+1, 0x1.8221d086914ebp+0, 0x1.670657e3a07bap-2,
     0x1.0f4951fd1e72dp-5, 0x1.b3874b8798286p-11};
  double cn0 = cn[0] + z2*cn[1];
  double cn2 = cn[2] + z2*cn[3];
  double cn4 = cn[4] + z2*cn[5];
  double cn6 = cn[6];
  cn0 += z4*cn2;
  cn4 += z4*cn6;
  cn0 += z8*cn4;
  cn0 *= z;
  double cd0 = cd[0] + z2*cd[1];
  double cd2 = cd[2] + z2*cd[3];
  double cd4 = cd[4] + z2*cd[5];
  double cd6 = cd[6];
  cd0 += z4*cd2;
  cd4 += z4*cd6;
  cd0 += z8*cd4;
  double r = cn0/cd0;
  if (gt) r = __builtin_copysign(pi2, z) - r;
  b64u64_u tr = {.f = r};
  u64 tail = (tr.u + 6)&(~0ul>>36);
  if(__builtin_expect(tail<=12, 0)){
    static const struct {union{float arg; unsigned uarg;}; float rh, rl;} st[] = {
      {{0x1.1ad646p-4f}, 0x1.1a6386p-4f, -0x1.fffffep-29f},
      {{0x1.f51a68p-11f}, 0x1.f51a5ep-11f, 0x1.ac7824p-62f},
      {{0x1.fc5d82p+0f}, 0x1.1ab2fp+0f, 0x1.0db9cap-52f},
      {{0x1.ddf9f6p+0f}, 0x1.143ec4p+0f, 0x1.5e8582p-54f},
      {{0x1.98c252p+12f}, 0x1.9215bp+0f, -0x1.069c58p-53f},
      {{0x1.71b3f4p+16f}, 0x1.921f04p+0f, -0x1.4d3ffcp-53f},
    };
    static const float q[] = {1.0f, -1.0f};
    for(int i=0;i<6;i++) {
      if(__builtin_expect(st[i].uarg == ax, 0))
	return  q[sgn]*st[i].rh + q[sgn]*st[i].rl;
    }
  }
  return r;
}

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

typedef union {float f; uint32_t u;} b32u32_u;
typedef union {double f; uint64_t u;} b64u64_u;
typedef uint64_t u64;

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
  /* now |x| >= 0x1p-13 */
  uint32_t ax = t.u&(~0u>>1);
  double z = x;
  if (gt) z = 1/z; /* gt is non-zero for |x| >= 1 */
  double z2 = z*z, z4 = z2*z2, z8 = z4*z4;
  /* polynomials generated using rminimax
     (https://gitlab.inria.fr/sfilip/rminimax) with the following command:
     ./ratapprox --function="atan(x)" --dom=[0.000122070,1] --num=[x,x^3,x^5,x^7,x^9,x^11,x^13] --den=[1,x^2,x^4,x^6,x^8,x^10,x^12] --output=atanf.sollya --log
     (see output atanf.sollya)
  */
  static const double cn[] =
    {0x1.51eccde075d67p-2, 0x1.a76bb5637f2f2p-1, 0x1.81e0eed20de88p-1,
     0x1.376c8ca67d11dp-2, 0x1.aec7b69202ac6p-5, 0x1.9561899acc73ep-9,
     0x1.bf9fa5b67e6p-16};
  static const double cd[] =
    {0x1.51eccde075d67p-2, 0x1.dfbdd7b392d28p-1, 0x1p+0,
     0x1.fd22bf0e89b54p-2, 0x1.d91ff8b576282p-4, 0x1.653ea99fc9bbp-7,
     0x1.1e7fcc202340ap-12};
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
    static const struct {union{float arg; uint32_t uarg;}; float rh, rl;} st[] = {
      {{0x1.1ad646p-4f}, 0x1.1a6386p-4f, -0x1.fffffep-29f},
      {{0x1.ddf9f6p+0f}, 0x1.143ec4p+0f, 0x1.5e8582p-54f},
      {{0x1.fc5d82p+0f}, 0x1.1ab2fp+0f, 0x1.0db9cap-52f},
    };
    static const float q[] = {1.0f, -1.0f};
    for(int i=0;i<3;i++) {
      if(__builtin_expect(st[i].uarg == ax, 0))
	return  q[sgn]*st[i].rh + q[sgn]*st[i].rl;
    }
  }
  return r;
}

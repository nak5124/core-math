/* Correctly-rounded tangent of binary32 value for angles in half-revolutions

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
#include <fenv.h>

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

/* __builtin_roundeven was introduced in gcc 10:
   https://gcc.gnu.org/gcc-10/changes.html,
   and in clang 17 */
#if (defined(__GNUC__) && __GNUC__ >= 10) || (defined(__clang__) && __clang_major__ >= 17)
#define HAS_BUILTIN_ROUNDEVEN
#endif

#ifndef HAS_BUILTIN_ROUNDEVEN
#include <math.h>
/* round x to nearest integer, breaking ties to even */
static float
__builtin_roundevenf (float x)
{
  float y = roundf (x); /* nearest, away from 0 */
  if (fabsf (y - x) == 0.5)
  {
    /* if y is odd, we should return y-1 if x>0, and y+1 if x<0 */
    union { float f; uint32_t n; } u, v;
    u.f = y;
    v.f = (x > 0) ? y - 1.0f : y + 1.0f;
    if (__builtin_ctz (v.n) > __builtin_ctz (u.n))
      y = v.f;
  }
  return y;
}
#endif

typedef union {float f; uint32_t u;} b32u32_u;

float cr_tanpif(float x){
  b32u32_u ix = {.f = x};
  uint32_t e = ix.u&(0xff<<23);
  if(__builtin_expect(e > (150<<23), 0)){
    if(e==(0xff<<23)){
      if(!(ix.u << 9)){
	errno = EDOM;
	feraiseexcept (FE_INVALID);
	return __builtin_nanf("inf");
      }
      return x + x; // nan
    }
    return __builtin_copysign(0.0f, x);
  }
  float x4 = 4.0f*x, nx4 = __builtin_roundevenf(x4), dx4 = x4-nx4;
  float ni = __builtin_roundevenf(x), zf = x-ni;
  if(__builtin_expect(dx4 == 0.0f, 0)){
    int k = x4;
    if(k&1) return __builtin_copysignf(1.0f,zf);
    k &= 7;
    if(k==0) return  __builtin_copysignf(0.0f,x);
    if(k==4) return  -__builtin_copysignf(0.0f,x);
    if(k==2) return  1.0f/0.0f;
    if(k==6) return -1.0f/0.0f;
  }
  ix.f = zf;
  uint32_t a = ix.u&(~0u>>1);
  if(__builtin_expect(a == 0x3e933802u, 0)) return  __builtin_copysignf(0x1.44cfbap+0f,zf) + __builtin_copysignf(0x1p-25f,zf);
  if(__builtin_expect(a == 0x38f26685u, 0)) return  __builtin_copysignf(0x1.7cc304p-12,zf) + __builtin_copysignf(0x1p-37f,zf);;
  double z = zf, z2 = z*z;

  static const double cn[] = {0x1.921fb54442d19p-1, -0x1.1f458b3e1f8d6p-2, 0x1.68a34bd0b8f6ap-6, -0x1.e4866f7a25f99p-13};
  static const double cd[] = {0x1p+0, -0x1.4b4b98d2df3a7p-1, 0x1.8e9926d2bb901p-4, -0x1.a6f77fd847eep-9};
  double z4 = z2*z2, r = (z-z*z2)*((cn[0]+z2*cn[1])+z4*(cn[2]+z2*cn[3]))/(((cd[0]+z2*cd[1])+z4*(cd[2]+z2*cd[3]))*(0.25-z2));
  return r;
}

#ifndef SKIP_C_FUNC_REDEF // icx provides this function
/* just to compile since glibc does not contain this function */
float tanpif(float x){
  return cr_tanpif(x);
}
#endif

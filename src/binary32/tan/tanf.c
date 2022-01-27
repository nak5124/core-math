/* Correctly-rounded tangent of binary32 value.

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
#include <fenv.h>
#include <errno.h>

/* __builtin_roundeven was introduced in gcc 10 */
#if defined(__GNUC__) && __GNUC__ >= 10
#define HAS_BUILTIN_ROUNDEVEN
#endif

#if defined(__clang__) && (defined(__AVX__) || defined(__SSE4_1__))
inline double __builtin_roundeven(double x){
   double ix;
#if defined __AVX__
   asm("vroundsd $0x8,%1,%1,%0":"=x"(ix):"x"(x));
#else /* __SSE4_1__ */
   asm("roundsd $0x8,%1,%0":"=x"(ix):"x"(x));
#endif
   return ix;
}
#define HAS_BUILTIN_ROUNDEVEN
#endif

#ifndef HAS_BUILTIN_ROUNDEVEN
#include <math.h>
/* round x to nearest integer, breaking ties to even */
static double
__builtin_roundeven (double x)
{
  double y = round (x); /* nearest, away from 0 */
  if (fabs (y - x) == 0.5)
  {
    /* if y is odd, we should return y-1 if x>0, and y+1 if x<0 */
    union { double f; uint64_t n; } u, v;
    u.f = y;
    v.f = (x > 0) ? y - 1.0 : y + 1.0;
    if (__builtin_ctz (v.n) > __builtin_ctz (u.n))
      y = v.f;
  }
  return y;
}
#endif

typedef union {float f; uint32_t u;} b32u32_u;
typedef union {double f; uint64_t u;} b64u64_u;
typedef unsigned __int128 u128;
typedef uint64_t u64;

static inline double rltl(float z, int *q){
  double x = z;
  double idl = -0x1.b1bbead603d8bp-32*x, idh = 0x1.45f306ep-1*x, id = __builtin_roundeven(idh);
  *q = (long)id;
  return (idh - id) + idl;
}

static double __attribute__((noinline)) rbig(uint32_t u, int *q){
  static const u64 ipi[] = {0xfe5163abdebbc562, 0xdb6295993c439041, 0xfc2757d1f534ddc0, 0xa2f9836e4e441529};
  int e = (u>>23)&0xff, i;
  u64 m = (u&(~0u>>9))|1<<23;
  u128 p0 = (u128)m*ipi[0];
  u128 p1 = (u128)m*ipi[1]; p1 += p0>>64;
  u128 p2 = (u128)m*ipi[2]; p2 += p1>>64;
  u128 p3 = (u128)m*ipi[3]; p3 += p2>>64;
  u64 p3h = p3>>64, p3l = p3, p2l = p2, p1l = p1;
  long a;
  int k = e-127, s = k-23;
  if(s<0){
    i =        p3h>>-s;
    a = p3h<<s|p3l>>-s;
  } else if(s==0) {
    i = p3h;
    a = p3l;
  } else if(s<64) {
    i = p3h<<s|p3l>>-s;
    a = p3l<<s|p2l>>-s;
  } else if(s==64) {
    i = p3l;
    a = p2l;
  } else {
    i = p3l<<s|p2l>>-s;
    a = p2l<<s|p1l>>-s;
  }
  int sgn = u; sgn >>= 31;
  long sm = a>>63;
  i -= sm;
  double z = (a^sgn)*0x1p-64;
  i = (i^sgn) - sgn;
  *q = i;
  return z;
}

float cr_tanf(float x){
  b32u32_u t = {.f = x};
  int e = (t.u>>23)&0xff, i;
  double z;
  if (__builtin_expect(e<127+28, 1)){
    if (__builtin_expect(e<115, 0)){
      if (__builtin_expect(e<102, 0))
	return __builtin_fmaf(x, __builtin_fabsf(x), x);
      float x2 = x*x;
      return __builtin_fmaf(x, 0x1.555556p-2f*x2, x);
    }
    z = rltl(x, &i);
  } else if (e<0xff){
    z = rbig(t.u, &i);
  } else {
    if(t.u<<9) return x; // nan
    errno = EDOM;
    feraiseexcept(FE_INVALID);
    return __builtin_nanf("tinf"); // inf
  }
  double z2 = z*z, z4 = z2*z2;
  static const double cn[] = {0x1.921fb54442d18p+0, -0x1.fd226e573289fp-2, 0x1.b7a60c8dac9f6p-6, -0x1.725beb40f33e5p-13};
  static const double cd[] = {0x1p+0, -0x1.2395347fb829dp+0, 0x1.2313660f29c36p-3, -0x1.9a707ab98d1c1p-9};
  static const double s[] = {0, 1};
  double n = cn[0] + z2*cn[1], n2 = cn[2] + z2*cn[3]; n += z4*n2;
  double d = cd[0] + z2*cd[1], d2 = cd[2] + z2*cd[3]; d += z4*d2;
  n *= z;
  double s0 = s[i&1], s1 = s[1-(i&1)];
  double r1 = (n*s1 - d*s0)/(n*s0 + d*s1);
  b64u64_u tr = {.f = r1};
  u64 tail = (tr.u + 7)&(~0ul>>35);
  if(__builtin_expect(tail<=14, 0)){
    static const struct {union{float arg; uint32_t uarg;}; float rh, rl;} st[] = {
      {{0x1.143ec4p+0f}, 0x1.ddf9f6p+0f, -0x1.891d24p-52f},
      {{0x1.ada6aap+27f}, 0x1.e80304p-3f, 0x1.419f46p-58f},
      {{0x1.af61dap+48f}, 0x1.60d1c8p-2f, -0x1.2d6c3ap-55f},
      {{0x1.0088bcp+52f}, 0x1.ca1edp+0f, 0x1.f6053p-53f},
      {{0x1.f90dfcp+72f}, 0x1.597f9cp-1f, 0x1.925978p-53f},
      {{0x1.cc4e22p+85f}, -0x1.f33584p+1f, 0x1.d7254ap-51f},
      {{0x1.a6ce12p+86f}, -0x1.c5612ep-1f, -0x1.26c33ep-53f},
      {{0x1.6a0b76p+102f}, -0x1.e42a1ep+0f, -0x1.1dc906p-52f},
    };
    uint32_t ax = t.u&(~0u>>1), sgn = t.u>>31;
    for(int i=0;i<8;i++) {
      if(__builtin_expect(st[i].uarg == ax, 0)){
	if(sgn)
	  return -st[i].rh - st[i].rl;
	else
	  return  st[i].rh + st[i].rl;
      }
    }
  }
  return r1;
}

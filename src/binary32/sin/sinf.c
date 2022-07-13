/* Correctly-rounded sine of binary32 value.

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
   __asm__("vroundsd $0x8,%1,%1,%0":"=x"(ix):"x"(x));
#else /* __SSE4_1__ */
   __asm__("roundsd $0x8,%1,%0":"=x"(ix):"x"(x));
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
    /* Negative shifts are undefined behaviour: p3l>>-s seems to work
       with gcc, but does not with clang. */
    i =        p3h>>(64-s);
    a = p3h<<s|p3l>>(64-s);
  } else if(s==0) {
    i = p3h;
    a = p3l;
  } else if(s<64) {
    i = p3h<<s|p3l>>(64-s);
    a = p3l<<s|p2l>>(64-s);
  } else if(s==64) {
    i = p3l;
    a = p2l;
  } else { /* s > 64 */
    i = p3l<<s|p2l>>(128-s);
    a = p2l<<s|p1l>>(128-s);
  }
  int sgn = u; sgn >>= 31;
  long sm = a>>63;
  i -= sm;
  double z = (a^sgn)*0x1p-64;
  i = (i^sgn) - sgn;
  *q = i;
  return z;
}

static inline double rltl(float z, int *q){
  double x = z;
  double idl = -0x1.b1bbead603d8bp-32*x, idh = 0x1.45f306ep-1*x, id = __builtin_roundeven(idh);
  *q = (long)id;
  return (idh - id) + idl;
}

float cr_sinf(float x){
  static const double
    cs[] = {
    -0x1.a51a6625307bdp-2, 0x1.9f9cb402b97ffp-5, -0x1.86a8e46de2fc4p-9,
    0x1.ac67ffda7836cp-14, -0x1.337d0b88f8cb7p-19, 0x1.3417d9749e139p-25},
    cc[] = {
      -0x1.3bd3cc9be458bp+0, 0x1.03c1f081b078ep-2, -0x1.55d3c7dbfe042p-6,
      0x1.e1f4fb610f151p-11, -0x1.a6c9c224d18abp-16, 0x1.f3dbf0909677fp-22},
    q[] = {1.0, 0, -1.0, 0};
  static const float tb[] = {1, -1};
  static const struct {union{float arg; uint32_t uarg;}; float rh, rl;} st[] = {
    {{0x1.fbd9c8p+22f}, -0x1.ff6dc2p-1f, 0x1.c23274p-57f},
    {{0x1.728fecp+37f}, -0x1.24f23cp-1f, 0x1.206be8p-54f},
    {{0x1.33333p+13f}, -0x1.63f4bap-2f,-0x1.fffffep-27f},
  };

  b32u32_u t = {.f = x};
  int e = (t.u>>23)&0xff, i;
  uint32_t ax = t.u&(~0u>>1), sgn = t.u>>31;
  double z;
  if (__builtin_expect(e<127+28, 1)){
    if (__builtin_expect(e<115, 0)){
      if (__builtin_expect(e<102, 0))
	return __builtin_fmaf(-x, __builtin_fabsf(x), x);
      float x2 = x*x;
      return __builtin_fmaf(-x, 0x1.555556p-3f*x2, x);
    }
    if(__builtin_expect(st[0].uarg == ax, 0)) return tb[sgn]*st[0].rh + tb[sgn]*st[0].rl;
    if(__builtin_expect(st[2].uarg == ax, 0)) return tb[sgn]*st[2].rh + tb[sgn]*st[2].rl;
    z = rltl(x, &i);
  } else if (e<0xff) {
    if(__builtin_expect(st[1].uarg == ax, 0)) return tb[sgn]*st[1].rh + tb[sgn]*st[1].rl;
    z = rbig(t.u, &i);
  } else {
    if(t.u<<9) return x; // nan
    errno = EDOM;
    feraiseexcept(FE_INVALID);
    return __builtin_nanf("cinf"); // inf
  }
  double z2 = z*z, z4 = z2*z2;
  double ms = q[i&3], mc = q[(i-1)&3];
  z *= 0x1.921fb54442d18p+0*ms;
  if(__builtin_expect(z2<0x1p-25, 0)){
    if(i&1){
      float a = z2;
      int j = (i>>1)&1;
      return __builtin_fmaf(((float)cc[0])*tb[j], a, tb[j]);
    }
  }
  double s0 = cs[0] + z2*cs[1];
  double s2 = cs[2] + z2*cs[3];
  double s4 = cs[4] + z2*cs[5];
  double rs = (z*z2)*(s0 + z4*(s2 + z4*s4)) + z;
  double c0 = cc[0] + z2*cc[1];
  double c2 = cc[2] + z2*cc[3];
  double c4 = cc[4] + z2*cc[5];
  double rc = (z2*mc)*(c0 + z4*(c2 + z4*c4)) + mc;
  return rs + rc;
}

/* Correctly-rounded atan2 function for two binary64 values.

Copyright (c) 2023 Paul Zimmermann

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

#include <stdio.h>
#include <stdint.h>

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union { double f; uint64_t u; } d64u64;
#include "tint.h"

// return non-zero if u (with sign bit cleared) encodes Inf or NaN
static inline int inf_or_nan (uint64_t u)
{
  return (u >> 52) == 0x7ff;
}

#define MASK 0x7ffffffffffffffful

static inline int is_nan (uint64_t u)
{
  u = u & MASK;
  uint64_t e = u >> 52;
  return e == 0x7ff && u != (e << 52);
}

static inline int is_inf (uint64_t u)
{
  return (u & MASK) == (0x7fful << 52);
}

// PI_H+PI_L approximates pi with error bounded by 2^-108.041
#define PI_H 0x1.921fb54442d18p+1
#define	PI_L 0x1.1a62633145c07p-53
// PI_OVER2_H+PI_OVER2_L approximates pi/2 with error bounded by 2^-109.041
#define PI_OVER2_H 0x1.921fb54442d18p+0
#define	PI_OVER2_L 0x1.1a62633145c07p-54
// PI_OVER4_H+PI_OVER4_L approximates pi/4 with error bounded by 2^-110.041
#define PI_OVER4_H 0x1.921fb54442d18p-1
#define	PI_OVER4_L 0x1.1a62633145c07p-55

/* The following is a degree-15 polynomial with odd coefficients
   approximating atan(x) on [0,2^-11.2] with maximal relative error
   2^-192.031 (cf atan2_small.sollya). */
static const tint_t Psmall[] = {
  // degree 1: 1
  {.h=0x8000000000000000, .m=0x0, .l=0x0, .ex=1, .sgn=0},
  // degree 3: -0x1.5555555555555555555555555555555555555555569de77ep-2
  {.h=0xaaaaaaaaaaaaaaaa, .m=0xaaaaaaaaaaaaaaaa, .l=0xaaaaaaaaab4ef3bf, .ex=-1, .sgn=1},
  // degree 5: 0x1.9999999999999999999999999999999999be8f1d48d2355cp-3
  {.h=0xcccccccccccccccc, .m=0xcccccccccccccccc, .l=0xccdf478ea4691aae, .ex=-2, .sgn=0},
  // degree 7: -0x1.2492492492492492492492492492c75ep-3
  {.h=0x9249249249249249, .m=0x24924924924963af, .l=0x0, .ex=-2, .sgn=1},
  // degree 9: 0x1.c71c71c71c71c71c71c71dbb1af83beap-4
  {.h=0xe38e38e38e38e38e, .m=0x38e38edd8d7c1df5, .l=0x0, .ex=-3, .sgn=0},
  // degree 11: -0x1.745d1745d1745d18p-4
  {.h=0xba2e8ba2e8ba2e8c, .m=0x0, .l=0x0, .ex=-3, .sgn=1},
  // degree 13: 0x1.3b13b13b13d919bap-4
  {.h=0x9d89d89d89ec8cdd, .m=0x0, .l=0x0, .ex=-3, .sgn=0},
  // degree 15: -0x1.1111103a0ee2178ep-4
  {.h=0x8888881d07710bc7, .m=0x0, .l=0x0, .ex=-3, .sgn=1},
};

static double
atan2_accurate_small (double y, double x)
{
  tint_t z[1], z2[1], p[1];
  printf ("y=%la z=%la\n", y, x);
  div_tint (z, y, x);
  printf ("z="); print_tint (z);
  mul_tint (z2, z, z);
  printf ("z2="); print_tint (z2);
  cp_tint (p, Psmall+7); // degree 15
  printf ("102: p="); print_tint (p);
  for (int i = 6; i >= 0; i--)
  {
    mul_tint (p, p, z2);
    printf ("106: p="); print_tint (p);
    add_tint (p, p, Psmall+i);
    printf ("108: p="); print_tint (p);
  }
  // multiply by z
  mul_tint (p, p, z);
  // printf ("p="); print_tint (p);
  return tint_tod (p);
}

// accurate path, assumes both y and x are neither NaN, nor +/-Inf, nor +/-0
static double
atan2_accurate (double y, double x)
{
  double z = y / x;
  if (__builtin_fabs (z) <= 0x1.bdb8cdadbe12p-12) // |z| < 2^-11.2
    return atan2_accurate_small (y, x);
  return 0;
}

// atan(y/x)
double cr_atan2 (double y, double x)
{
  d64u64 uy = {.f = y}, ux = {.f = x};
  uint64_t ay = uy.u & MASK, ax = ux.u & MASK;

  if (__builtin_expect (inf_or_nan (ay) || inf_or_nan (ax), 0))
  {
    if (is_nan (ay) || is_nan (ax))
      return y + x; // if y or x is sNaN, returns qNaN are raises invalid
    // Now neither y nor x is NaN, but at least one is +Inf or -Inf
    if (is_inf (ay) && is_inf (ax)) // both y and x are +/-Inf
    {
      // atan2 (+/-Inf,-Inf) = +/-3pi/4
      if (x < 0)
        return (y > 0) ? 3 * PI_OVER4_H + 3 * PI_OVER4_L
          : -3 * PI_OVER4_H - 3 * PI_OVER4_L;
      // atan2 (+/-Inf,+Inf) = +/-pi/4
      return (y > 0) ? PI_OVER4_H + PI_OVER4_L : -PI_OVER4_H - PI_OVER4_L;
    }
    // now only one of y and x is +/-Inf
    if (is_inf (ax))
    {
      if (x < 0)
        return (uy.u >> 63) ? -PI_H - PI_L : PI_H + PI_L;
      // atan2(+/-0,x) = +/-0 for x > 0
      // atan2(+/-y,+Inf) = +/-0 for finite y>0
      return __builtin_copysign (0, y);
    }
    // now y = +/-Inf
    // atan2(+/-Inf,x) = +/-pi/2 for finite x
    return (y > 0) ? PI_OVER2_H + PI_OVER2_L : -PI_OVER2_H - PI_OVER2_L;
  }

  if (__builtin_expect (y == 0 || x == 0, 0))
  {
    if (y == 0 && x == 0)
    {
      if (ux.u == 0) // atan2(+/-0, +0) = +/-0
        return y;
      // atan2(+/-0, +0) = +/-pi
      return (uy.u == 0) ? PI_H + PI_L : -PI_H - PI_L;
    }
    // only one of y and x is zero
    if (y == 0)
    {
      // atan2(+/-0,x) = +/-0 for x>0
      if (x > 0) return y;
      // atan2(+/-0,x) = +/-pi for x<0
      return (uy.u == 0) ? PI_H + PI_L : -PI_H - PI_L;
    }
    // now only x is zero
    // atan2(y,+/-0) = -pi/2 for y<0
    // atan2(y,+/-0) = +pi/2 for y>0
    return (y > 0) ? PI_OVER2_H + PI_OVER2_L : -PI_OVER2_H - PI_OVER2_L;
  }

  // now both y and x are neither NaN, nor +/-Inf, nor +/-0

  return atan2_accurate (y, x);

  return 0.0;
}

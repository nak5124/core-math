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

// PI_OVER2_H+PI_OVER2_L approximates pi/2 with error bounded by 2^-109.041
#define PI_OVER2_H 0x1.921fb54442d18p+0
#define	PI_OVER2_L 0x1.1a62633145c07p-54

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
        return (y > 0) ? 3 * PI_OVER2_H / 2 + 3 * PI_OVER2_L / 2
          : -3 * PI_OVER2_H / 2 - 3 * PI_OVER2_L / 2;
      // atan2 (+/-Inf,+Inf) = +/-pi/4
      return (y > 0) ? PI_OVER2_H / 2 + PI_OVER2_L / 2
        : -PI_OVER2_H / 2 - PI_OVER2_L / 2;
    }
    // now only one of y and x is +/-Inf
    if (is_inf (ax))
    {
      if (x < 0)
        return (uy.u >> 63) ? -2 * PI_OVER2_H - 2 * PI_OVER2_L
          : +2 * PI_OVER2_H + 2 * PI_OVER2_L;
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
      return (uy.u == 0) ? 2 * PI_OVER2_H + 2 * PI_OVER2_L
        : -2 * PI_OVER2_H - 2 * PI_OVER2_L;
    }
    // only one of y and x is zero
    if (y == 0)
    {
      // atan2(+/-0,x) = +/-0 for x>0
      if (x > 0) return y;
      // atan2(+/-0,x) = +/-pi for x<0
      return (uy.u == 0) ? 2 * PI_OVER2_H + 2 * PI_OVER2_L
        : -2 * PI_OVER2_H - 2 * PI_OVER2_L;
    }
    // now only x is zero
    // atan2(y,+/-0) = -pi/2 for y<0
    // atan2(y,+/-0) = +pi/2 for y>0
    return (y > 0) ? PI_OVER2_H + PI_OVER2_L : -PI_OVER2_H - PI_OVER2_L;
  }

  return 0.0;
}

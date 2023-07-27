/* Correctly-rounded expm1 function for binary64 value.

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

#include <stdint.h>

typedef union {double f; uint64_t u;} b64u64_u;

/* Given -0x1.2b708872320e2p+5 < x < -0x1.6a09e667f3bccp-53 or
   0x1.6a09e667f3bccp-53 < x < 0x1.62e42fefa39fp+9, put in h + l a
   double-double approximation of expm1(x), and return the maximal
   corresponding absolute error. */
static double
expm1_fast (double *h, double *l, double x)
{
  *h = 0;
  *l = 0;
  return 1;
}

double
cr_expm1 (double x)
{
  b64u64_u t = {.f = x};
  uint64_t ux = t.u, ax = ux & 0x7ffffffffffffffflu;

  if (__builtin_expect (ux >= 0xc042b708872320e2, 0))
  {
    // x = -NaN or x <= -0x1.2b708872320e2p+5
    if ((ux >> 52) == 0xfff) // -NaN
      return x;
    // for x <= -0x1.2b708872320e2p+5, expm1(x) rounds to -1 to nearest
    return -1.0 + 0x1p-53;
  }
  else if (__builtin_expect (ax >= 0x40862e42fefa39f0, 0))
  {
    // x = +NaN or x >= 0x1.62e42fefa39fp+9
    if ((ux >> 52) == 0x7ff) // +NaN
      return x;
    // for x >= 0x1.62e42fefa39fp+9, expm1(x) rounds to +Inf to nearest
    return 0x1.fffffffffffffp+1023 + 0x1.fffffffffffffp+1023;
  }
  else if (ax <= 0x3ca6a09e667f3bcc) // |x| <= 0x1.6a09e667f3bccp-53
    /* then expm1(x) rounds to x (to nearest), with Taylor expansion
       x + x^2/2 + ... */
  {
    if (ax < 0x3ca0000000000000)
    {
      /* |x| < 0x1p-53: x^2 < 1/2 ulp(x), we have to deal with -0 apart
         since fma (-0, -0, -0) is (+0) + (-0) which evaluates to +0
         for some rounding modes */
      return (x == 0) ? x : __builtin_fma (x, x, x);
    }
    else
      /* 0x1p-53 <= |x| <= 0x1.6a09e667f3bccp-53: x/4 is exactly
         representable, and x^2/4 < 1/2 ulp(x) */
      return __builtin_fma (x, x * 0.25, x);
  }

  /* -0x1.2b708872320e2p+5 < x < -0x1.6a09e667f3bccp-53 or
     0x1.6a09e667f3bccp-53 < x < 0x1.62e42fefa39fp+9 */

  double err, h, l;
  err = expm1_fast (&h, &l, x);
  double left = h + (l - err), right = h + (l + err);
  if (left == right)
    return left;

  return 0;
}

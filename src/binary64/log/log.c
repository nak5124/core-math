/* Correctly rounded logarithm of binary64 values.

Copyright (c) 2022 Paul Zimmermann, Inria.

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
#include <math.h>

typedef union { double f; uint64_t u; } d64u64;

/* Add a + b exactly, such that *hi + *lo = a + b.
   Assumes |a| >= |b|.  */
static void
fast_two_sum (double *hi, double *lo, double a, double b)
{
  double e;

  *hi = a + b;
  e = *hi - a; /* exact */
  *lo = b - e; /* exact */
  /* Now hi + lo = a + b exactly for rounding to nearest.
     For directed rounding modes, this is not always true.
     Take for example a = 1, b = 2^-200, and rounding up,
     then hi = 1 + 2^-52, e = 2^-52 (it can be proven that
     e is always exact), and lo = -2^52 + 2^-105, thus
     hi + lo = 1 + 2^-105 <> a + b = 1 + 2^-200.
     A bound on the error is given
     in "Tight interval inclusions with compensated algorithms"
     by Stef Graillat and Fabienne Jézéquel,
     IEEE Transactions on Computers, 2020. Proposition 3.2 says that
     the difference between a+b and hi+lo is bounded by 4u^2|a+b|
     and also by 4u^2|hi|. Here u=2^-53, thus we get:
     |(a+b)-(hi+lo)| <= 2^-104 min(|a+b|,|hi|) */
}

/* given 1 <= x < 2, put in h+l a double-double approximation of log(x) */
static void
cr_log_fast (double *h, double *l, double x)
{
  d64u64 v = {.f = x};
  int i = (v >> 44) & 0xff; /* 0 <= i < 256 */
}

double
cr_log (double x)
{
  if (x <= 0.0)
  {
    /* f(x<0) is NaN, f(+/-0) is -Inf and raises DivByZero */
    if (x < 0)
      return 0.0 / 0.0;
    else
      return 1.0 / -0.0;
  }
  /* now x > 0 */
  d64u64 v = {.f = x};
  int e = (v.u >> 52) - 0x3ff, bias = 0;
  if (e == 0x400) /* +Inf or NaN */
    return x;
  /* now 0 < x < +Inf */
  if (e == -0x3ff)
  {
    v.f *= 0x1p52;
    bias = 52;
    e = (v.u >> 52) - 0x3ff;
  }
  v.u -= (int64_t) e << 52;
  double m = v.f;
  e -= bias;
  /* now x = m*2^e with 1 <= m < 2 */
  double h, l;
  cr_log_fast (&h, &l, m);
  if (e != 0)
  {
    /* Add e*log(2) to (h,l), where -1074 <= e <= 1023, thus e has at most
       11 bits. We store log2_h on 42 bits, so that e*log2_h is exact. */
    static double log2_h = 0x1.62e42fefa38p-1, log2_l = 0x1.ef35793c7673p-45;
    /* |log(2) - (h+l)| < 2^-102.01 */
    double hh = e * log2_h; /* exact */
    double ll = __builtin_fma (e, log2_l, l);
    fast_two_sum (&h, &l, hh, h); /* rounding error bounded by 2^-104*|hh| */
    l += ll;
  }
  /* Error bound:
     a) error from cr_log_fast: XXX
     b) approx. error of log(2) multiplied by e: 1074*2^-102.01 < 2^-91.94
     c) rounding error on e*log2_l+l: assuming |ll| <= 2^-33, we get
        ulp(2^-34) = 2^-86
     Thus we get 2^-85.97 + XXX.
   */
  return h + l;
}

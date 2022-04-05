/* Correctly rounded exponential function for binary64 values.

Copyright (c) 2021-2022 Paul Zimmermann, Inria.

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

#include <math.h>
#include <fenv.h>
#include <assert.h>

/* Add a + b exactly, such that *hi + *lo = a + b.
   Assumes |a| >= |b| and rounding to nearest.  */
static void
fast_two_sum (double *hi, double *lo, double a, double b)
{
  double e;

  *hi = a + b;
  e = *hi - a; /* exact */
  *lo = b - e; /* exact */
  /* Now *hi + *lo = a + b exactly.  */
}

/* h + l <- a * b */
static void
dekker (double *h, double *l, double a, double b)
{
  *h = a * b;
#ifdef __FP_FAST_FMA
  *l = __builtin_fma (a, b, -*h);
#else /* use Dekker's algorithm */
#define MAGIC 0x8000001
  double ah = a * MAGIC, bh = b * MAGIC;
  ah = (a - ah) + ah;
  bh = (b - bh) + bh;
  double al = a - ah;
  double bl = b - bh;
  *l = (((ah * bh - *h) + ah * bl) + al * bh) + al * bl;
#endif
}

double
cr_exp (double x)
{
  if (isnan (x))
    return x + x; /* always return qNaN, even for sNaN input */

  /* Other special cases, cf page 1390, left column.  */

  /* Deal with overflow.  */
#define X_MAX 0x1.62e42fefa39efp+9
  if (x > X_MAX)
  {
    /* The +0x1p970 returns 0x1.fffffffffffffp1023 or +Inf according to the
       rounding mode.  */
    double x0 = 0x1.fffffffffffffp1023; /* largest representable number */
    /* If x > x0, then x is +Inf.  */
    return (x > x0) ? x : x0 + 0x1p970;
  }
  
  /* deal with underflow */
#define X_MIN -0x1.6232bdd7abcd2p+9
#define X_DMIN -0x1.74385446d71c3p+9
#define X_DMIN2 -0x1.74910d52d3051p+9
  if (x < X_MIN)
  {
    if (x < X_DMIN && rnd != FE_TONEAREST)
    {
      /* Warning: 0x1p-1074 * 0.5 is optimized by gcc to 0, whatever the
         rounding mode: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=34678.
         We add 'volatile' which seems to solve this issue. */
      volatile double xmin = 0x1p-1074;
      return xmin * 0.5;
    }
    if (x < X_DMIN2 && rnd == FE_TONEAREST)
      return 0;
  }

  /* now -0x1.6232bdd7abcd2p+9 <= x <= 0x1.62e42fefa39efp+9 */
  
  /* first multiply x by a double-double approximation of 1/log(2) */
  static const double log2_h = 0x1.71547652b82fep+0;
  static const double log2_l = 0x1.777d0ffda0ep-56;
  double h, l;
  dekker (&h, &l, x, log2_h);
  l += x * log2_l;
  
  /* now x/log(2) ~ h + l thus exp(x) ~ 2^h * 2^l where |l| < 2^-42 */
}

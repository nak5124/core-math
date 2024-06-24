/* Correctly rounded reciprocal square root for long-double values.

Copyright (c) 2024 Alexei Sibidanov and Paul Zimmermann

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

#define TRACE 0x1.8003bb2d200784p-16385L

#include <stdint.h>
#include <fenv.h>

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

// anonymous structs, see https://port70.net/~nsz/c/c11/n1570.html#6.7.2.1p19
typedef union {
  long double f;
  struct __attribute__((__packed__))
  {uint64_t m; uint32_t e:16; uint32_t empty:16;};
} b96u96_u;

typedef union {double f;uint64_t u;} b64u64_u;

// Multiply exactly a and b, such that *hi + *lo = a * b.
static inline void a_mul_double (double *hi, double *lo, double a, double b) {
  *hi = a * b;
  *lo = __builtin_fma (a, b, -*hi);
}

/* put in (h+l)*2^e an approximation of 1/sqrt(x) where x = vm/2^63*2^(e-16383)
   with |h + l - 1/sqrt(xr))| < 2^-97.654, where xr is the reduced operand,
   1/2 <= xr < 2. */
static void
fast_path (double *h, double *l, unsigned long vm, int *e)
{
  /* convert vm/2^63 to a double-double representation xh + xl */
  b64u64_u th = {.u = (0x3fful<<52) | (vm >> 11)};
  b64u64_u tl = {.u = (0x3cbul<<52) | ((vm << 53) >> 12)};
  double xh = th.f, xl = tl.f - 0x1p-52;
  // 1 <= xh < 2 and 0 <= xl < 2^-52

  *e -= 16383; // unbias e

  // if e is odd, divide xh,xl by 2
  if (*e & 1)
  {
    xh = xh * 0.5;
    xl = xl * 0.5;
    (*e) ++;
  }
  
  *e = - (*e / 2);
  
  // 1/sqrt(x) = 1/sqrt(xh+xl)*2^e with 1/2 <= xh, xh + xl < 2

  double yh, yl;
  yh = 1.0 / __builtin_sqrt (xh); // 1/2 <= yh <= 2
  /* Let s == __builtin_sqrt (xh), we have s = sqrt(xh) * (1 + eps1)
     with |eps1| < 2^-52.
     Then yh = 1/s * (1 + eps2) with |eps2| < 2^-52, thus
          yh = 1/sqrt(xh) * (1 + eps2)/(1 + eps1)
             = 1/sqrt(xh) * (1 + eps3) with |eps3| < 2^-50.999 */

  /* Perform one step of Newton iteration: y' = y - y/2 * (x * y^2 - 1).
     Let e = x*y^2-1 and e' = x*y'^2-1. Since y' = y - y/2*e, we deduce:
     y'^2 = y^2 - y^2*e + y^2/4*e^2
     x*y'^2-1 = (x*y^2-1) - x*y^2*e + x*y^2/4*e^2
     e' = e - x*y^2*e + x*y^2/4*e^2
        = (1-x*y^2)*e + x*y^2/4*e^2
        = e^2 + x*y^2/4*e^2
        = e^2 * (1 + (e+1)/4)
     Using y=yh, we obtain:
     e = (xh+xl)/xh * (1+eps3)^2 - 1
       = (1+eps3)^2 - 1 + xl/xh*(1+eps3)^2
     thus since |xl/xh| < 2^-52:
     |e| <= 2^-49.677.
     If we plug in e' = e^2 * (1 + (e+1)/4) we obtain:
     |e'| < 2^-99.032. */

  double zh, zl;
  a_mul_double (&zh, &zl, yh, yh); // zh+zl = yh^2, exact
  // since 1/2 <= yh <= 2, 1/4 <= zh+zl <= 4, with |zl| < ulp(zh) <= 2^-51
  // x * y^2 - 1 = (x * zh - 1) + x * zl
  yl = __builtin_fma (xh, zh, -1.0);
  /* since yh = 1/sqrt(xh) * (1 + eps3), yh^2 = 1/xh * (1+eps3)^2,
     thus (zh+zl)*xh = (1+eps3)^2
     |zh*xh-1| <= |zl*xh| + (1+eps3)^2 - 1
               <= 2^-51*2 + 2^-49.998 <= 2^-48.998
     we deduce that |yl| < 2^-48, and the rounding error of the above fma
     is bounded by ulp(2^-48.998) = 2^-101. */
  yl = __builtin_fma (xh, zl, yl);
  /* |xh| <= 2, |zl| <= 2^-51, and |yl| <= 2^-48.998 thus the new value of
     yl is < 2*2^-51 + 2^-48.998 <= 2^-48.413, and the rounding error of the
     above fma is bounded by ulp(2^-48.413) = 2^-101 again. */
  // add xl * zh
  yl = __builtin_fma (xl, zh, yl);
  /* |xl| <= 2^-52, |zh| <= 4 and |yl| <= 2^-48.413, thus the new value of
     yl is < 2^-52*4 + 2^-48.413 <= 2^-47.998, and the rounding error of the
     above fma is bounded by ulp(2^-47.998) = 2^-100.
     We neglected the term xl*zl which is bounded by 2^-52*2^-51, thus
     induces an error of at most 2^-103.
     The absolute error on yl is bounded by:
     2^-101 + 2^-101 + 2^-100 + 2^-103 < 2^-98.912. */
  yl = yh * yl * -0.5;
  /* since |yh| <= 2 and |yl| <= 2^-47.998, the new value of yl is bounded by
     2^-47.998 too, and the rounding error in this last multiplication (the
     multiplication by -0.5 is exact) is bounded by ulp(2^-47.998) = 2^-100.
     The total error is bounded by:
     * 2^-99.032 for the mathematical error e' (see above)
     * 2^-98.912 for the error induced by the rounding error on the previous
       value of yl (multiplied by yh with |yh| <= 2 and by 0.5)
     * 2^-100 for the rounding error in this last operation
     This gives an absolute error bounded by:

     |yh + yl - 1/sqrt(xh + xl)| < 2^-99.032 + 2^-98.912 + 2^-100 < 2^-97.654
  */

  *h = yh;
  *l = yl;
}

long double
cr_rsqrtl (long double x)
{
  b96u96_u v = {.f = x};
  int e = v.e & 0x7fff;

  // if (x == TRACE) printf ("cr_rsqrtl: x=%La e=%d\n", x, e);

  // check NaN, Inf, 0
  if (__builtin_expect (x < 0 || e == 32767 || (e == 0 && v.m == 0), 0))
  {
    if (x == 0) return 1.0L / x;   // x=+0 and x=-0
    if (x < 0) return 0.0L / 0.0L; // x<0: rsqrt(x)=NaN
    if (x > 0) return +0L;         // x=Inf
    return x;                      // x=NaN
  }

  //  if (x == TRACE) printf ("v.m=%lx e=%d\n", v.m, e);

  // rsqrt(x) is exact iff x = 2^(2k)
  if (__builtin_expect (e == 0 || (v.m == 0x8000000000000000 && (e & 1)), 0))
  {
    if (e > 0) // normal numbers
    {
      // for x=1, e=16383
      v.e = 16383 + (16383 - e) / 2;
      return v.f;
    }
    // case e=0: subnormal numbers
    // x = 2^(2k) iff v.m == 2^(2t+1)
    int cnt = __builtin_ctzll (v.m);
    if ((cnt & 1) && (v.m == ((uint64_t) 1 << cnt)))
    {
      v.m = 0x8000000000000000ul;
      // x = 2^(-16445+cnt)
      v.e = 16383 + (16445 - cnt) / 2;
      return v.f;
    }
    // normalize subnormal numbers not of the form 2^(2k)
    cnt = __builtin_clzll (v.m);
    v.m <<= cnt;
    e -= cnt - 1;
  }

  double h, l;
  fast_path (&h, &l, v.m, &e);
  // if (x == TRACE) printf ("h=%la l=%la\n", h, l);
  long double H = h, L = l;
  const long double err = 0x1.46p-98L; // 2^-97.654 < 0x1.46p-98
  long double left = H + (L - err), right = H + (L + err);
  // printf ("left=%La right=%La\n", left, right);
  if (__builtin_expect (left == right, 1))
    return __builtin_ldexpl (left, e);

  return x;
}

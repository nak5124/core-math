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

#define TRACE 0x0.000000000000002p-16385L

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
   with |h + l - 1/sqrt(x))| < XXX */
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
  
  // 1/sqrt(x) = 1/sqrt(xh+xl)*2^e with 1/2 <= xh + xl < 2

  double yh, yl;
  yh = 1.0 / __builtin_sqrt (xh);

  // perform one step of Newton iteration: y' = y - y/2 * (x * y^2 - 1)
  double zh, zl;
  a_mul_double (&zh, &zl, yh, yh); // zh+zl = yh^2
  // x * y^2 - 1 = (x * zh - 1) + x * zl
  yl = __builtin_fma (xh, zh, -1.0);
  yl = __builtin_fma (xl, zh, yl);
  yl = __builtin_fma (xh, zl, yl);
  yl = yh * yl * -0.5;

  // yh+yl approximates 1/sqrt(xh+xl) with about 106 bits of accuracy
  *h = yh;
  *l = yl;
}

long double
cr_rsqrtl (long double x)
{
  b96u96_u v = {.f = x};
  int e = v.e & 0x7fff;

  // check NaN, Inf, 0
  if (__builtin_expect (e == 32767 || (e == 0 && v.m == 0), 0))
  {
    if (x == 0) return 1.0L / x; // +0 and -0
    if (x < 0) return x / x;     // -Inf
    if (x > 0) return +0L;       // Inf
    return x;                    // NaN
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
  }

  double h, l;
  fast_path (&h, &l, v.m, &e);
  long double H = h, L = l;
  const long double err = 0x1p-100L;
  long double left = H + (L - err), right = H + (L + err);
  // printf ("left=%La right=%La\n", left, right);
  if (__builtin_expect (left == right, 1))
    return __builtin_ldexpl (left, e);

  return x;
}

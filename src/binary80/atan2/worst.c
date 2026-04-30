/* Generate worst cases for atan2l

Copyright (c) 2023-2026 Paul Zimmermann, Inria.

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

Usage:

   $ gcc -O3 worst.c -o worst -lmpfr -lm -fopenmp
   $ ./worst <nnn>

   where <nnn> is the wanted minimal number of identical bits after the round
   bit.
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <mpfr.h>
#include <math.h>
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif
#include <fenv.h>
#include <unistd.h>

double maxerr;
unsigned long nsols = 0;

typedef union { __int128 n; long double x; } union_t;

static inline __int128
asuint80 (long double x)
{
  union_t u;
  u.x = x;
  return u.n;
}

static void
check_one (__int128 p, __int128 q, mpfr_exp_t e)
{
  mpfr_t x, y, z, t;
  mpfr_init2 (x, 64);
  mpfr_init2 (y, 64);
  mpfr_init2 (z, 65);
  mpfr_init2 (t, 200);
  mpfr_set_ui_2exp (y, (uint64_t) p, e, MPFR_RNDN);
  if (q >= 0)
    mpfr_set_ui (x, q, MPFR_RNDN);
  else {
    mpfr_set_ui (x, -q, MPFR_RNDN);
    mpfr_neg (x, x, MPFR_RNDN);
  }
  mpfr_atan2 (z, y, x, MPFR_RNDN);
  double err;
  while (1)
  {
    mpfr_atan2 (t, y, x, MPFR_RNDN);
    mpfr_sub (t, t, z, MPFR_RNDN);
    if (mpfr_cmp_ui (t, 0) != 0)
    {
      mpfr_exp_t e = mpfr_get_exp (z);
      // ulp(z) = 2^(e-65)
      mpfr_mul_2exp (t, t, 65 - e, MPFR_RNDN);
      err = fabs (mpfr_get_d (t, MPFR_RNDN));
      break;
    }
    mpfr_set_prec (t, 2 * mpfr_get_prec (t));
  }
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp critical
#endif
  if (err < maxerr)
  {
    mpfr_printf ("%Ra,%Ra\n", y, x);
    fflush (stdout);
    nsols ++;
    if (nsols == 1000) exit (0);
  }
  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
  mpfr_clear (t);
}

/* given z, find the largest convergent p/q of tan(z)/2^e such that
   p,q < 2^64, with 1/2 <= |tan(z)|/2^e < 1, and output the corresponding
   worst case if it has at least m identical bits after the round bit.
   If add_one_ulp is non-zero, add 1/2 ulp to z (above). */
static void
confrac (long double z, int add_one_ulp)
{
  mpfr_prec_t prec;
  mpfr_t yl, yh;
  int sign, e, ok = 0, last_i = -1;

  prec = 200;
  while (ok == 0) {
  ok = 1;
  mpfr_init2 (yl, 65);
  mpfr_init2 (yh, prec);
  mpfr_set_ld (yl, z, MPFR_RNDN);
  if (add_one_ulp)
    mpfr_nextabove (yl);
  mpfr_prec_round (yl, prec, MPFR_RNDN);
  mpfr_tan (yl, yl, MPFR_RNDD);
  sign = (mpfr_sgn (yl) > 0) ? 1 : -1;
  mpfr_abs (yl, yl, MPFR_RNDN);
  e = mpfr_get_exp (yl); /* 2^(e-1) <= yl < 2^e */
  mpfr_div_2si (yl, yl, e, MPFR_RNDN);
  /* now 1/2 <= yl < 1 */
  mpfr_set (yh, yl, MPFR_RNDN);
  mpfr_nextabove (yh);
  /* yl <= |tan(z)|/2^e < yh */
  uint64_t p0 = 1, p1 = 0, q0 = 0, q1 = 1, p2, q2;
  for (int i = 0; ; i++)
  {
    mpfr_ui_div (yl, 1, yl, MPFR_RNDU);
    mpfr_ui_div	(yh, 1, yh, MPFR_RNDD);
    mpfr_swap (yl, yh);
    /* yl <= 1/y <= yh */
    uint64_t al = mpfr_get_ui (yl, MPFR_RNDD);
    uint64_t ah = mpfr_get_ui (yh, MPFR_RNDD);
    if (al != ah)
    {
      ok = 0;
      prec = 2 * prec;
      goto end;
    }
    mpfr_sub_ui (yl, yl, al, MPFR_RNDD);
    mpfr_sub_ui (yh, yh, ah, MPFR_RNDU);
    p2 = al * p1 + p0;
    q2 = al * q1 + q0;
    __int128 p2i = (__int128) al * (__int128) p1 + (__int128) p0;
    __int128 q2i = (__int128) al * (__int128) q1 + (__int128) q0;
    /* in case prec was increased, we don't want to check twice the
       same convergents */
    if (i > last_i)
    {
      if ((p2i >> 64) || (q2i >> 64))
      {
        check_one (p1, (sign == 1) ? (__int128) q1 : - (__int128) q1, e);
        break;
      }
      last_i = i;
    }
    p0 = p1; p1 = p2;
    q0 = q1; q1 = q2;
  }
  end:
  mpfr_clear (yl);
  mpfr_clear (yh);
  }
}

static void
doit (long double z)
{
  confrac (z, 0);
  confrac (z, 1);
}

int
main (int argc, char *argv[])
{
  assert (argc == 2);
  int m = atoi (argv[1]);
  maxerr = ldexp (1.0, -m);
  //  srand48 (getpid ());
  while (1)
  {
#define LENGTH 1024
    long double buf[LENGTH];
    for (int i = 0; i < LENGTH; i++)
    {
      long double d = drand48 ();
      long double e = drand48 ();
      buf[i] = 0.5l + d * 0.5l + e * 0x1p-54l;
    }
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for schedule(static,1)
#endif
      for (int i = 0; i < LENGTH; i++)
        doit (buf[i]);
  }
}

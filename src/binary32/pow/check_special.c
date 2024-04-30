/* Generate special cases for powf testing.

Copyright (c) 2022-2023 Stéphane Glondu and Paul Zimmermann, Inria.

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
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#include <mpfr.h>
#include <math.h>
#include <omp.h>
#include <assert.h>

float cr_powf (float, float);
float ref_pow (float, float);
int ref_fesetround (int);
void ref_init (void);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

int mid = 1; // if mid=1, also check midpoint cases

static inline uint32_t
asuint (float f)
{
  union
  {
    float f;
    uint32_t i;
  } u = {f};
  return u.i;
}

static void
check (float x, float y)
{
  float z1, z2;
  ref_init();
  ref_fesetround(rnd);
  mpfr_flags_clear (MPFR_FLAGS_INEXACT);
  z1 = ref_pow(x, y);
  mpfr_flags_t inex1 = mpfr_flags_test (MPFR_FLAGS_INEXACT);
  fesetround(rnd1[rnd]);
  feclearexcept (FE_INEXACT);
  z2 = cr_powf(x, y);
  fexcept_t inex2;
  fegetexceptflag (&inex2, FE_INEXACT);
  if (asuint (z1) != asuint (z2)) {
    printf("FAIL x=%a y=%a ref=%a z=%a\n", x, y, z1, z2);
    fflush(stdout);
    exit(1);
  }
  if ((inex1 == 0) && (inex2 != 0))
  {
    printf ("Spurious inexact exception for x=%a y=%a (z=%a)\n", x, y, z1);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
  if ((inex1 != 0) && (inex2 == 0))
  {
    printf ("Missing inexact exception for x=%a y=%a (z=%a)\n", x, y, z1);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
}

/* y positive integer: both x and -x are solutions */
static void
count_uint_y (int y)
{
  float xmin, xmax;
  int emin, emax;
  mpfr_t z;
  mpfr_init2 (z, 25);
  mpfr_set_ui_2exp (z, 1, -150, MPFR_RNDN);
  mpfr_rootn_ui (z, z, y, MPFR_RNDU);
  xmin = mpfr_get_flt (z, MPFR_RNDU);
  mpfr_set_ui_2exp (z, 1, 128, MPFR_RNDN);
  mpfr_nextbelow (z);
  mpfr_rootn_ui (z, z, y, MPFR_RNDD);
  xmax = mpfr_get_flt (z, MPFR_RNDD);
  /* we want xmin <= x <= xmax, with x^y exact */
  /* write x = m*2^e with m odd, m >= 3 */
  /* compute the maximum odd m such that m^y < 2^(24+mid) */
  mpfr_set_ui_2exp (z, 1, 24 + mid, MPFR_RNDN);
  mpfr_nextbelow (z);
  mpfr_rootn_ui (z, z, y, MPFR_RNDD);
  assert (mpfr_fits_sint_p (z, MPFR_RNDD));
  int maxm = mpfr_get_ui (z, MPFR_RNDD);
  /* since m*2^e <= xmax, we have 2^e <= xmax/m <= xmax/3 */
  mpfr_set_flt (z, xmax, MPFR_RNDD);
  mpfr_div_ui (z, z, 3, MPFR_RNDD);
  emax = mpfr_get_exp (z) - 1;
  /* for the minimal exponent, since x is an odd multiple of 2^e,
     x^y is an odd multiple of 2^(y*e), thus we want y*e >= -149-mid */
  emin = - ((149 + mid) / y);
  for (int e = emin; e <= emax; e++)
    for (int m = 3; m <= maxm; m += 2)
    {
      mpfr_set_ui_2exp (z, m, e, MPFR_RNDN);
      if (mpfr_cmp_d (z, xmin) >= 0 && mpfr_cmp_d (z, xmax) <= 0)
      {
        float x = mpfr_get_flt (z, MPFR_RNDN);
        check (x, (float) y);
        check (-x, (float) y);
      }

    }
  mpfr_clear (z);
}

/* y = n/2^f */
static void
count_uint_2exp_y (int n, int f)
{
  int F = 1 << f;
  float xmin, xmax;
  int emin, emax;
  mpfr_t z;
  float y = ldexpf (n, -f);
  mpfr_init2 (z, 25);
  mpfr_set_ui_2exp (z, 1, -150, MPFR_RNDN);
  mpfr_rootn_ui (z, z, n, MPFR_RNDU);
  mpfr_pow_ui (z, z, F, MPFR_RNDU);
  xmin = mpfr_get_flt (z, MPFR_RNDU);
  mpfr_set_ui_2exp (z, 1, 128, MPFR_RNDN);
  mpfr_nextbelow (z);
  mpfr_rootn_ui (z, z, n, MPFR_RNDD);
  mpfr_pow_ui (z, z, F, MPFR_RNDD);
  xmax = mpfr_get_flt (z, MPFR_RNDD);
  /* we want xmin <= x <= xmax, with x^y exact */
  /* write x = m*2^e with m odd, m >= 3, m = k^(2^f) */
  /* we should have both k^(2^f) < 2^24 and k^n < 2^(24+mid) */
  /* compute the maximum odd k such that k^(2^f) < 2^24 */
  mpfr_set_ui_2exp (z, 1, 24, MPFR_RNDN);
  mpfr_nextbelow (z);
  mpfr_rootn_ui (z, z, F, MPFR_RNDD);
  assert (mpfr_fits_sint_p (z, MPFR_RNDD));
  int maxk = mpfr_get_ui (z, MPFR_RNDD);
  /* compute the maximum odd k such that k^n < 2^(24+mid) */
  mpfr_set_ui_2exp (z, 1, 24 + mid, MPFR_RNDN);
  mpfr_nextbelow (z);
  mpfr_rootn_ui (z, z, n, MPFR_RNDD);
  assert (mpfr_fits_sint_p (z, MPFR_RNDD));
  int maxk2 = mpfr_get_ui (z, MPFR_RNDD);
  maxk = (maxk < maxk2) ? maxk : maxk2;
  /* Write x=m*2^e, we should have m=k^(2^f) and e multiple of 2^f,
     then x^y = k^n*2^(e*n/2^f).
     We should have m=k^(2^f) with k <= kmax. */
  /* since m*2^e <= xmax, we have 2^e <= xmax/m <= xmax/3 */
  mpfr_set_flt (z, xmax, MPFR_RNDD);
  mpfr_div_ui (z, z, 3, MPFR_RNDD);
  emax = mpfr_get_exp (z) - 1;
  /* for the minimal exponent, since x is an odd multiple of 2^e,
     x^y is an odd multiple of 2^(y*e), thus we want e >= -149 and
     y*e >= -149-mid */
  emin = - ((149 + mid) / y);
  if (emin < -149)
    emin = -149; /* so that x is representable */
  /* we should have e multiple of F */
  while (emin % F)
    emin ++;
  for (int e = emin; e <= emax; e += F)
  {
    for (int k = 3; k <= maxk; k += 2)
    {
      unsigned long m = k;
      for (int j = 0; j < f; j++)
        m = m * m;
      /* m (odd) should be less than 2^24 */
      assert (m < 0x1000000);
      mpfr_set_ui_2exp (z, m, e, MPFR_RNDN);
      if (mpfr_cmp_d (z, xmin) >= 0 && mpfr_cmp_d (z, xmax) <= 0)
        check (mpfr_get_flt (z, MPFR_RNDN), y);
    }
  }
  mpfr_clear (z);
}

/* x = +/-2^e:
   for y integer, and -149-mid <= e*y <= 127, both x=2^e and -2^e are solutions
   for y=n/2^f for f>=1, n odd, we need e divisible by 2^f,
   and -149-mid <= e*y/2^f <= 127 */
static void
count_pow2_x (int e)
{
  if (e == 0) /* trivial solutions */
    return;

  /* case y integer */
  float x = ldexpf (1.0f, e);
  int ymin, ymax;
  if (e > 0)
  {
    ymin = - ((149+mid) / e);
    ymax = 127 / e;
  }
  else /* e < 0 */
  {
    ymin = - (127 / -e);
    ymax = (149+mid) / (-e);
  }
  for (int y = ymin; y <= ymax; y++)
    if (y != 0 && y != 1)
    {
      check (x, (float) y);
      check (-x, (float) y);
    }

  /* case y = n/2^f */
  int f = 1;
  while ((e % 2) == 0)
  {
    /* invariant: e = e_orig/2^f */
    e = e / 2;
    /* -149-mid <= e*y <= 127 */
    if (e > 0)
    {
      ymin = - ((149+mid) / e);
      ymax = 127 / e;
    }
    else /* e < 0 */
    {
      ymin = - (127 / -e);
      ymax = (149+mid) / (-e);
    }
    /* y should be odd */
    if ((ymin % 2) == 0)
      ymin ++;
    for (int y = ymin; y <= ymax; y += 2)
      check (x, ldexpf (y, -f));
    f ++;
  }
}

// check exact and midpoint cases
static void
check_exact (void)
{
  /* First deal with integer y >= 2. If x is not a power of 2, then y <= 15
     whatever the value of mid, since 3^15 has 24 bits, and 3^16 has 26 bits.
     Indeed, assume x = m*2^e with m odd, then m >= 3, thus we should have
     m^y < 2^(24+mid). */
  for (int n = 2; n <= 15; n++)
    count_uint_y (n);
  /* Now deal with y = n/2^f for a positive integer f, and an odd n. If x is
     not a power of 2, assume x = m*2^e with m odd, m >= 3. Then m^(1/2^f)
     should be an odd integer >= 3, which implies f <= 3 whatever the value of
     mid, since 3^(2^3) has 13 bits, and 3^(2^4) has 26 bits.
     For the same reason as above, n <= 15. */
  for (int f = 1; f <= 3; f++)
    for (int n = 1; n <= 15; n += 2)
      count_uint_2exp_y (n, f);
  /* Now deal with x=2^e. */
  for (int e = -149; e <= 127; e++)
    count_pow2_x (e);
}

int
main (int argc, char *argv[])
{
  while (argc >= 2)
    {
      if (strcmp (argv[1], "--rndn") == 0)
        {
          rnd = 0;
          argc --;
          argv ++;
        }
      else if (strcmp (argv[1], "--rndz") == 0)
        {
          rnd = 1;
          argc --;
          argv ++;
        }
      else if (strcmp (argv[1], "--rndu") == 0)
        {
          rnd = 2;
          argc --;
          argv ++;
        }
      else if (strcmp (argv[1], "--rndd") == 0)
        {
          rnd = 3;
          argc --;
          argv ++;
        }
      else if (strcmp (argv[1], "--verbose") == 0)
        {
          verbose = 1;
          argc --;
          argv ++;
        }
      else
        {
          fprintf (stderr, "Error, unknown option %s\n", argv[1]);
          exit (1);
        }
    }

  printf ("Checking exact and midpoint values\n");
  fflush (stdout);
  check_exact ();
  return 0;
}

/* Generate special cases for atan2pif testing.

Copyright (c) 2022 Stéphane Glondu and Paul Zimmermann, Inria.

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
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <fenv.h>
#include <math.h>
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif
#include <mpfr.h>

float cr_atan2pif (float, float);
void ref_init (void);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };
int rnd2[] = { MPFR_RNDN, MPFR_RNDZ, MPFR_RNDU, MPFR_RNDD };

int rnd = 0;

int verbose = 0;

/* reference code using MPFR */
static float
ref_atan2pi (float y, float x)
{
  mpfr_t xi, yi;
  mpfr_inits2 (24, xi, yi, NULL);
  mpfr_set_flt (xi, x, MPFR_RNDN);
  mpfr_set_flt (yi, y, MPFR_RNDN);
  int inex = mpfr_atan2pi (xi, yi, xi, rnd2[rnd]);
  mpfr_subnormalize (xi, inex, rnd2[rnd]);
  float ret = mpfr_get_flt (xi, MPFR_RNDN);
  mpfr_clears (xi, yi, NULL);
  return ret;
}

typedef union { uint32_t n; float x; } union_t;

static float
asfloat (uint32_t n)
{
  union_t u;
  u.n = n;
  return u.x;
}

static void
check (float y, float x)
{
  float z, t;
  mpfr_flags_clear (MPFR_FLAGS_INEXACT);
  t = ref_atan2pi (y, x);
  mpfr_flags_t inex1 = mpfr_flags_test (MPFR_FLAGS_INEXACT);
  feclearexcept (FE_INEXACT);
  z = cr_atan2pif (y, x);
  fexcept_t inex2;
  fegetexceptflag (&inex2, FE_INEXACT);
  if ((isnan (t) && !isnan(z)) || (!isnan (t) && isnan(z)) ||
      (!isnan (t) && !isnan(z) && z != t))
  {
    printf ("FAIL y=%a x=%a ref=%a z=%a\n", y, x, t, z);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
#ifdef CORE_MATH_CHECK_INEXACT
  if ((inex1 == 0) && (inex2 != 0))
  {
    printf ("Spurious inexact exception for x=%a y=%a\n", x, y);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
  if ((inex1 != 0) && (inex2 == 0))
  {
    printf ("Missing inexact exception for x=%a y=%a\n", x, y);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
#endif
}

#define N 100000000

static void
check_random (int i)
{
  long l;
  float x, y;
  struct drand48_data buffer[1];
  ref_init ();
  fesetround (rnd1[rnd]);
  srand48_r (i, buffer);
  for (int n = 0; n < N; n++)
  {
    lrand48_r (buffer, &l);
    y = asfloat (l);
    lrand48_r (buffer, &l);
    x = asfloat (l);
    check (y, x);
    check (y, -x);
    check (-y, x);
    check (-y, -x);
  }
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

  int nthreads = 1;
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel
  nthreads = omp_get_num_threads ();
#endif
  /* check random values */
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (int i = 0; i < nthreads; i++)
    check_random (getpid () + i);
  return 0;
}

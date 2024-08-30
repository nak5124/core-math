/* Check correctness of binary32 function like sincos by exhaustive search.

Copyright (c) 2022 Alexei Sibidanov.
Copyright (c) 2022-2024 Paul Zimmermann, INRIA.

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
#include <math.h>
#include <mpfr.h>
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif

#include "function_under_test.h"

void cr_function_under_test (float, float*, float*);
void ref_function_under_test (float, float*, float*);
int ref_fesetround (int);
void ref_init (void);

/* the code below is to check correctness by exhaustive search */

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int keep = 0;

typedef union { uint32_t n; float x; } union_t;

float
asfloat (uint32_t n)
{
  union_t u;
  u.n = n;
  return u.x;
}

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

static int
is_equal (float y1, float y2)
{
  if (isnan (y1))
    return isnan (y2);
  if (isnan (y2))
    return isnan (y1);
  return asuint (y1) == asuint (y2);
}

void
doit (uint32_t n)
{
  float x, y1, y2, z1, z2;
  x = asfloat (n);
  ref_init ();
  ref_fesetround (rnd);
  mpfr_flags_clear (MPFR_FLAGS_INEXACT);
  ref_function_under_test (x, &y1, &y2);
  mpfr_flags_t inex_y = mpfr_flags_test (MPFR_FLAGS_INEXACT);
  fesetround (rnd1[rnd]);
  feclearexcept (FE_INEXACT);
  cr_function_under_test (x, &z1, &z2);
  fexcept_t inex_z;
  fegetexceptflag (&inex_z, FE_INEXACT);
  if (!is_equal (y1, z1) || !is_equal (y2, z2))
  {
    printf ("FAIL x=%a ref=(%a,%a) z=(%a,%a)\n", x, y1, y2, z1, z2);
    fflush (stdout);
    if (!keep) exit (1);
  }
#ifdef CORE_MATH_CHECK_INEXACT
  if ((inex_y == 0) && (inex_z != 0))
  {
    printf ("Spurious inexact exception for x=%a\n", x);
    fflush (stdout);
    if (!keep) exit (1);
  }
  if ((inex_y != 0) && (inex_z == 0))
  {
    printf ("Missing inexact exception for x=%a\n", x);
    fflush (stdout);
    if (!keep) exit (1);
  }
#endif
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
      else if (strcmp (argv[1], "--keep") == 0)
        {
          keep = 1;
          argc --;
          argv ++;
        }
      else
        {
          fprintf (stderr, "Error, unknown option %s\n", argv[1]);
          exit (1);
        }
    }

  return doloop();
}

/* Generate special cases for expl testing.

Copyright (c) 2022-2024 Stéphane Glondu and Paul Zimmermann, Inria.

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
#include <unistd.h>
#include <math.h> // for ldexpl

#include "cm_types.h"

int ref_fesetround (int);
void ref_init (void);

long double cr_expl (long double);
long double ref_expl (long double);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

static int
is_nan (long double x)
{
  b80u80_t v = {.f = x};
  return ((v.e == 0x7fff || v.e == 0xffff) && (v.m != (1ul << 63)));
}

static inline int
is_equal (long double x, long double y)
{
  if (is_nan (x))
    return is_nan (y);
  if (is_nan (y))
    return is_nan (x);
  return x == y;
}

static void
check (long double x)
{
  long double y1 = ref_expl (x);
  fesetround (rnd1[rnd]);
  long double y2 = cr_expl (x);
  if (! is_equal (y1, y2))
  {
    printf ("FAIL x=%La ref=%La z=%La\n", x, y1, y2);
    fflush (stdout);
    exit (1);
  }
}

static long double
get_random ()
{
  b80u80_t v;
  v.m = rand ();
  v.m |= (uint64_t) rand () << 31;
  v.m |= (uint64_t) (rand () & 1) << 62;
	// the low 63 bits of m are random
  v.e = rand () & 0xffff;
  // if e is not 0 nor 0x8000 (0 or subnormal), m should have its most
  // significant bit set, otherwise it should be cleared
  // cf https://en.wikipedia.org/wiki/Extended_precision
  uint64_t t = (v.e&0x7fff) != 0;
  v.m |= t << 63;
  return v.f;
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

  ref_init();
  ref_fesetround (rnd);

  unsigned int seed = getpid ();
  srand (seed);

#define N 10000000UL /* total number of tests */

  printf ("Checking results in subnormal range\n");
  /* check subnormal results */
  /* x0 is the smallest x such that 2^-16446 <= RN(exp(x)) */
  long double x0 = -0x1.643bfcfe13c57552p+13L;
  /* x1 is the smallest x such that 2^-16384 <= RN(exp(x) */
  long double x1 = -0x1.62e42fefa39ef356p+13L;
  /* in the [x0,x1) range, floating-point numbers have ulp = 2^-50 */
  long double ulp = 0x1p-50L; // ulp = ulp(x0) = ulp(x1)
  long double dx = (x1 - x0) / (long double) N; // total numbers in [x0,x1]
  unsigned long skip = dx / ulp; // distance between two checked numbers
  int n0 = seed % skip;
  x0 += (long double) n0 * ulp;  // we start at a random x0
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (uint64_t n = 0; n < N; n++)
    check (x0 + (long double) n * dx);

  /* x2 is the smallest x such that 2^-16382 <= RN(exp(x)) */
  long double x2 = -0x1.62d918ce2421d65ep+13L;
  dx = (x2 - x1) / (long double) N;
  x1 += (long double) n0 * ulp;  // we start at a random x1
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (uint64_t n = 0; n < N; n++) {
    ref_init ();
    ref_fesetround (rnd);
    check (x1 + (long double) n * dx);
	}

	printf("Checking results near overflow\n");
	/*x3 is the biggest x such that exp(x) < MAX_LDBL*/
	long double x3 = 0x1.62e42fefa39ef357p+13L;
	dx             = 0x0.0000000010000000p+13L/ (long double) N;
	ulp = 0x1p-50L;
	skip = dx / ulp;
	n0 = seed % skip;
	x3 += (long double) n0 * ulp - (skip/2) * ulp;
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
	for (uint64_t n = 0; n < N; n++) {
    ref_init ();
    ref_fesetround (rnd);
		check(x3 + (long double) n * dx);
	}

  printf ("Checking random values with |x| < 2^-20\n");
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (uint64_t n = 0; n < N; n++)
  {
    ref_init ();
    ref_fesetround (rnd);
    long double x;
    x = get_random ();
    int e;
    x = frexpl (x, &e);
    check (ldexpl (x, -20));
  }

  printf ("Checking random values\n");
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (uint64_t n = 0; n < N; n++)
  {
    ref_init ();
    ref_fesetround (rnd);
    long double x;
    x = get_random ();
    check (x);
  }

  return 0;
}

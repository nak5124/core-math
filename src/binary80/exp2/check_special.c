/* Generate special cases for exp2l testing.

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
#include <math.h>
#include <unistd.h>

int ref_fesetround (int);
void ref_init (void);

double cr_exp2l (double);
double ref_exp2l (double);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

static inline int
is_equal (long double x, long double y)
{
  if (isnan (x))
    return isnan (y);
  if (isnan (y))
    return isnan (x);
  return x == y;
}

static void
check (long double x)
{
  double y1 = ref_exp2l (x);
  fesetround (rnd1[rnd]);
  double y2 = cr_exp2l (x);
  if (! is_equal (y1, y2))
  {
    printf ("FAIL x=%La ref=%La z=%La\n", x, y1, y2);
    fflush (stdout);
    exit (1);
  }
}

typedef union {long double f; uint64_t u, v;} b80u80_u;

static long double
get_random ()
{
  b80u80_u v;
  v.u = rand ();
  v.u |= (uint64_t) rand () << 31;
  v.u |= (uint64_t) rand () << 62;
  v.v = rand ();
  v.v |= (uint64_t) rand () << 31;
  v.v |= (uint64_t) rand () << 62;
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

  printf ("Checking results in subnormal range\n");
  /* check subnormal results */
  /* x0 is the smallest x such that 2^-16446 <= RN(2^-x) */
  long double x0 = 16446;
  /* x1 is the smallest x such that 2^-16384 <= RN(2^-x) */
  double x1 = 16384;
  /* in the [x0,x1) range, floating-point numbers have an integer part
     of 15 bits, thus we multiply by 2^49 to get integers, where 49 = 64-15 */
  uint64_t n0 = ldexpl (x0, 49);
  uint64_t n1 = ldexpl (x1, 49);
#define SKIP 320000
  n0 += getpid () % SKIP;
#pragma omp parallel for
  for (uint64_t n = n0; n < n1; n += SKIP)
    check (-ldexpl ((long double) n, -49));
  /* x2 is the smallest x such that 2^-16382 <= RN(2^-x) */
  double x2 = 16382;
  /* in the [x1,x2) range, floating-point numbers have an integer part
     of 14 bits, thus we multiply by 2^50 to get integers, where 50 = 64-14 */
  n1 = ldexpl (x1, 50);
  uint64_t n2 = ldexpl (x2, 50);
#pragma omp parallel for
  for (uint64_t n = n1; n < n2; n += SKIP)
    check (-ldexpl ((lonog double) n, -50));

  printf ("Checking random values\n");
#define N 1000000000UL /* total number of tests */

  unsigned int seed = getpid ();
  srand (seed);

#pragma omp parallel for
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

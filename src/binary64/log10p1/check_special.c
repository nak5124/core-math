/* Check log10p1 on random inputs.

Copyright (c) 2022-2023 Paul Zimmermann, Inria.

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
#include <sys/types.h>
#include <unistd.h>
#include <omp.h>

int ref_init (void);
int ref_fesetround (int);

double cr_log10p1 (double);
double ref_log10p1 (double);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

static inline uint64_t
asuint64 (double f)
{
  union
  {
    double f;
    uint64_t i;
  } u = {f};
  return u.i;
}

typedef union {double f; uint64_t u;} b64u64_u;

static double
get_random ()
{
  b64u64_u v;
  long l;
  v.u = rand ();
  v.u |= (uint64_t) rand () << 31;
  v.u |= (uint64_t) rand () << 62;
  return v.f;
}

static void
check (double x)
{
  int bug;
  double y1 = ref_log10p1 (x);
  fesetround (rnd1[rnd]);
  double y2 = cr_log10p1 (x);
  if (isnan (y1))
    bug = !isnan (y2);
  else if (isnan (y2))
    bug = !isnan (y1);
  else
    bug = asuint64 (y1) != asuint64 (y2);
  if (bug)
  {
    printf ("FAIL x=%la ref=%la z=%la\n", x, y1, y2);
    fflush (stdout);
    exit (1);
  }
}

/* check the inexact exception is not raised when log10p1(x) is exact */
static void
check_inexact (void)
{
  double x = 0, y;
  fexcept_t flagp;
  /* log10p1(x) is exact for x = 10^i-1 for 0 <= i <= 15
     (10^16-1 has 54 bits and is odd) */
  for (int i = 0; i <= 15; i++)
  {
    feclearexcept (FE_INEXACT);
    y = cr_log10p1 (x);
    fegetexceptflag (&flagp, FE_INEXACT);
    if (flagp)
    {
      printf ("Inexact flag set for x=%la\n", x);
      fflush (stdout);
      exit (1);
    }
    x = 10.0 * x + 9.0;
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
  ref_init ();
  ref_fesetround (rnd);

  printf ("Checking inexact flag\n");
  check_inexact ();

#define N 1000000000UL /* total number of tests */

  unsigned int seed = getpid ();
  srand (seed);

  printf ("Checking random values\n");
#pragma omp parallel for
  for (uint64_t n = 0; n < N; n++)
  {
    ref_init ();
    ref_fesetround (rnd);
    double x;
    do x = get_random (); while (x <= -1.0);
    check (x);
  }

  return 0;
}

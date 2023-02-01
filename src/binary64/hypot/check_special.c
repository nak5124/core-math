/* Generate special cases for hypot testing.

Copyright (c) 2022-2023 Stéphane Glondu, Paul Zimmermann, Inria.

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
#include <omp.h>
#include <unistd.h>

void doloop (int, int);
extern double cr_hypot (double, double);
extern int ref_fesetround (int);
extern void ref_init (void);
extern mpfr_rnd_t rnd2[];
extern double ref_hypot (double, double);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd;
int verbose = 0;

typedef union {double f; uint64_t u;} b64u64_u;

static double
get_random (struct drand48_data *buffer)
{
  b64u64_u v;
  long l;
  lrand48_r (buffer, &l);
  v.u = l;
  lrand48_r (buffer, &l);
  v.u |= (uint64_t) l << 31;
  lrand48_r (buffer, &l);
  v.u |= (uint64_t) l << 62;
  return v.f;
}

static void
check (double x, double y)
{
  double z, t;
  mpfr_t X, Y, Z;
  mpfr_init2 (X, 53);
  mpfr_init2 (Y, 53);
  mpfr_init2 (Z, 53);
  mpfr_set_d (X, x, MPFR_RNDN);
  mpfr_set_d (Y, y, MPFR_RNDN);
  z = cr_hypot (x, y);
  t = ref_hypot (x, y);
  if (z != t)
  {
    printf ("cr_hypot and ref_hypot differ for x=%la y=%la\n", x, y);
    printf ("cr_hypot  gives %la\n", z);
    printf ("ref_hypot gives %la\n", t);
    exit (1);
  }
  mpfr_clear (X);
  mpfr_clear (Y);
  mpfr_clear (Z);
}

static void
check_random (int i)
{
  ref_init ();
  ref_fesetround (rnd);
  fesetround(rnd1[rnd]);
  struct drand48_data buffer[1];
  double x, y;
  srand48_r (i, buffer);
  while (1)
  {
    x = get_random (buffer);
    y = get_random (buffer);
    x = 0x1.282a03d9ba3bap-727;
    y = 0x1.cf90629a60c16p-711;
    check (x, y);
  }
}

static void
check_random_all (void)
{
  int nthreads;
#pragma omp parallel
  nthreads = omp_get_num_threads ();
#pragma omp parallel for
  for (int i = 0; i < nthreads; i++)
    check_random (getpid () + i);
}

int
main (int argc, char *argv[])
{
  int random = 0;
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
      else if (strcmp (argv[1], "--random") == 0)
        {
          random = 1;
          argc --;
          argv ++;
        }
      else
        {
          fprintf (stderr, "Error, unknown option %s\n", argv[1]);
          exit (1);
        }
    }

  if (random)
    check_random_all ();

  /* we check triples with exponent difference 0 <= k <= 26 */
  doloop(0, 26);
  return 0;
}

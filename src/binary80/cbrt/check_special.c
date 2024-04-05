/* Generate special cases for cbrtl testing.

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
#include <math.h> // for ldexpl, cbrtl

int ref_fesetround (int);
void ref_init (void);

long double cr_cbrtl (long double);
long double ref_cbrtl (long double);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

// only the lower 16 bits of e are used
// 1.0 has encoding m=2^63, e=16383
// -1.0 has encoding m=2^63, e=49151
// 2 has encoding m=2^63, e=16384
// +qnan has encoding m=2^63+2^62, e=32767
// -qnan has encoding m=2^63+2^62, e=65535
// +inf has encoding m=2^63, e=32767
// -inf has encoding m=2^63, e=65535
// +snan has encoding m=2^63+2^62-1, e=32767
// -snan has encoding m=2^63+2^62-1, e=65535
typedef union {long double f; struct {uint64_t m; uint16_t e;};} b80u80_t;

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
  b80u80_t v = {.f = x}, w = {.f = y};
  return v.e == w.e && v.m == w.m; // ensures +0 and -0 differ
}

static void
check (long double x)
{
  long double y1 = ref_cbrtl (x);
  fesetround (rnd1[rnd]);
  long double y2 = cr_cbrtl (x);
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
  v.m |= (uint64_t) rand () << 62;
  v.e = rand () & 0xffff;
  // if e is not 0 nor 0x7fff nor 0xffff, m should have its msb set
  uint64_t t = v.e != 0 && v.e != 0x7fff && v.e != 0xffff;
  v.m |= t << 63;
  return v.f;
}

// check exact values (m*2^e)^3 with 2^61 <= |m^3| < 2^64
static void
check_exact (void)
{
  // the smallest exact cube is 2^-16443, it is generated with m=2097152 and e=-5502;
  // the largest exact cube is 0x1.ffffdbd247267c7ap+16383, it is generated with
  // m=2642245 and e=5440
#pragma omp parallel for
  for (int e = -5502; e <= 5440; e++)
  {
    uint64_t m;
    // 2^61 <= |m^3| < 2^64 implies 1321123 <= m <= 2642245
    for (m = 1321123; m <= 2642245; m++)
    {
      long double x;
      x = ldexpl ((long double) m, e);
      x = x * x * x;
      check (x);
      check (-x);
    }
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

  ref_init();
  ref_fesetround (rnd);

  printf ("Checking exact values\n");
  check_exact ();

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

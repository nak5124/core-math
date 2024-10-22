/* Generate special cases for hypotf testing.

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
#include <unistd.h>
#include <fenv.h>
#include <math.h>
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif
#include <assert.h>

#include "cm_types.h"

int ref_fesetround (int);
void ref_init (void);

// triples.c
extern void doloop (int, int);
extern void check (float, float);
extern uint64_t gcd (uint64_t, uint64_t);
float cr_hypotf (float, float);

// worst_p1.c
extern void doit_subnormal_above (uint32_t);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

/* Check all Pythagorean triples z^2 = x^2 + y^2 with z in the subnormal
   range. We necessarily have x = r^2 - s^2, y = 2*r*s, z = r^2 + s^2
   with gcd(r,s) = 1 and one of r, s even
   (see https://oeis.org/wiki/Pythagorean_triples).
*/
static void
check_triples_subnormal (void)
{
  uint64_t r, s, x, y, z;
  /* the smallest denormal is 2^-149, the smallest normal is 2^-126,
     thus x, y, z are of the form k*2^-149 with k < 2^23. */

  // type I: r is odd
  for (r = 1; r <= 2896; r += 2)
    for (s = 2; s < r; s += 2)
    {
      if (gcd (r, s) == 1)
      {
        x = r * r - s * s;
        y = 2 * r * s;
        z = r * r + s * s;
        if (z > 0x7fffff)
          break;
        // now (x,y,z) is a primitive Pythagorean triple
        for (int n = 1; ; n++)
        {
          uint64_t nn = n * n;
          uint64_t xx = x * nn, yy = y * nn, zz = z * nn;
          if (zz > 0x7fffff)
            break;
          check (ldexpf (xx, -149), ldexpf (yy, -149));
        }
      }
    }

  // type II: r is even
  for (r = 2; r <= 2896; r += 2)
    for (s = 1; s < r; s += 2)
    {
      if (gcd (r, s) == 1)
      {
        x = r * r - s * s;
        y = 2 * r * s;
        z = r * r + s * s;
        if (z > 0x7fffff)
          break;
        // now (x,y,z) is a primitive Pythagorean triple
        for (int n = 1; ; n++)
        {
          uint64_t nn = n * n;
          uint64_t xx = x * nn, yy = y * nn, zz = z * nn;
          if (zz > 0x7fffff)
            break;
          check (ldexpf (xx, -149), ldexpf (yy, -149));
        }
      }
    }
}

/* Check pairs (x,y) in subnormal range such that x = u*2^-149, y = v*2^-149,
   with u^2 + v^2 = w^2 + 1, u <= v. We force 2 <= u to avoid the trivial
   solutions u=1, v=w. See https://oeis.org/A050796. */
static void
check_triples_subnormal_above (void)
{
  doit_subnormal_above (8388608);
}

static float
get_random (struct drand48_data *buffer)
{
  b32u32_u v;
  int64_t l;
  lrand48_r (buffer, &l);
  v.u = l;
  // lrand48_r generates only 31 bits
  lrand48_r (buffer, &l);
  v.u |= (uint32_t) l << 31;
  return v.f;
}

#define N 10000000ul

static void
check_random (int i)
{
  ref_init ();
  ref_fesetround (rnd);
  fesetround(rnd1[rnd]);
  struct drand48_data buffer[1];
  float x, y;
  srand48_r (i, buffer);
  for (uint64_t n = 0; n < N; n++)
  {
    x = get_random (buffer);
    y = get_random (buffer);
    check (x, y);
  }
}

static void
check_random_all (void)
{
  int nthreads = 1;
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel
  nthreads = omp_get_num_threads ();
#endif
#pragma omp parallel for
  for (int i = 0; i < nthreads; i++)
    check_random (getpid () + i);
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

  printf ("Checking random values\n");
  check_random_all ();

  /* we check triples with exponent difference 0 <= k <= 12 */
  printf ("Checking near-exact subnormal values\n");
  fflush (stdout);
  check_triples_subnormal_above ();
  printf ("Checking exact subnormal values\n");
  fflush (stdout);
  check_triples_subnormal ();
  printf ("Checking Pythagorean triples\n");
  fflush (stdout);
  doloop(0, 12);
  return 0;
}

/* Generate special cases for log2l testing.

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
#include <mpfr.h>
#include <assert.h>
#include <unistd.h>
#include <math.h> // for ldexpl

#include "cm_types.h"

int ref_fesetround (int);
void ref_init (void);

long double cr_log2l (long double);
long double ref_log2l (long double);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

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
  uint64_t t = (v.e & 0x7fff) != 0;
  v.m |= t << 63;
  return v.f;
}

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
  mpfr_flags_clear (MPFR_FLAGS_INEXACT);
  long double y1 = ref_log2l (x);
  mpfr_flags_t inex1 = mpfr_flags_test (MPFR_FLAGS_INEXACT);
  fesetround (rnd1[rnd]);
  feclearexcept (FE_INEXACT);
  long double y2 = cr_log2l (x);
  fexcept_t inex2;
  fegetexceptflag (&inex2, FE_INEXACT);
  if (!is_equal (y1, y2))
  {
    printf ("FAIL x=%La ref=%La z=%La\n", x, y1, y2);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }
#ifdef CORE_MATH_CHECK_INEXACT
  if ((inex1 == 0) && (inex2 != 0))
  {
    printf ("Spurious inexact exception for x=%La (y=%La)\n", x, y1);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }
  if ((inex1 != 0) && (inex2 == 0))
  {
    printf ("Missing inexact exception for x=%La (y=%La)\n", x, y1);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
#endif
}

/* check all x*2^k such that log2(x*2^k) and log2(x) are in the same binade,
   i.e., k+log2(x) and log2(x) are in the same binade */
static void
check_extended (long double x)
{
  if (is_nan (x) || 2 * x == x) // avoid NaN, +/-Inf, and +/-0
    return;

  assert (x > 0);

  int e;
  long double r = frexpl (x, &e);
  // x = r*2^e with 0.5 <= r < 1
  
  /* For x < 1, log2(r*2^e) and log2(r*2^f) are in the same binade
     if -2^(k+1) < e, f <= -2^k, since -1 <= log2(r) < 0 has the same sign
     as e, f.
     For x > 1, log2(r*2^e) and log2(r*2^f) are in the same binade
     if 2^k < e, f <= 2^(k+1), since -1 <= log2(r) < 0 has sign opposite
     to that of e, f. */

  int f = (e >= 0) ? e : -e;
  int K = 1;
  while (f >= 2 * K)
    K = 2 * K;
  // now K <= f < 2*K [except for e=0]
  int emin, emax;
  if (e == 0 || e == 1)
    emin = emax = e;
  else if (e < 0) {
    K = -K;
    // e <= -1 thus K <= -1
    emin = 2 * K + 1;
    emax = K;
  } else {
    // e >= 2 thus K >= 2
    // 2^k = e is not allowed
    if (K == e)
      K = K / 2;
    emin = K + 1;
    emax = 2 * K;
  }
  assert (emin <= e && e <= emax);
  // if emin < -16445, then r*2^e is rounded to zero (to nearest)
  if (emin < -16445)
    emin = -16445;
  for (e = emin; e <= emax; e++)
    check (ldexpl (r, e));
}

static void
readstdin(long double **result, int *count)
{
  char *buf = NULL;
  size_t buflength = 0;
  ssize_t n;
  int allocated = 512;

  *count = 0;
  if (NULL == (*result = malloc(allocated * sizeof(long double)))) {
    fprintf(stderr, "malloc failed\n");
    exit(1);
  }

  while ((n = getline(&buf, &buflength, stdin)) >= 0) {
    if (n > 0 && buf[0] == '#') continue;
    if (*count >= allocated) {
      int newsize = 2 * allocated;
      long double *newresult = realloc(*result, newsize * sizeof(long double));
      if (NULL == newresult) {
        fprintf(stderr, "realloc(%d) failed\n", newsize);
        exit(1);
      }
      allocated = newsize;
      *result = newresult;
    }
    long double *item = *result + *count;
    // special code for snan, since glibc does not read them
    if (strncmp (buf, "snan", 4) == 0 || strncmp (buf, "+snan", 5) == 0)
    {
      b80u80_t v;
      // +snan has encoding m=2^63+1, e=32767 (for example)
      v.e = 0x7fff;
      v.m = 0x8000000000000001ul;
      *item = v.f;
      (*count)++;
    }
    else if (strncmp (buf, "-snan", 5) == 0)
    {
      b80u80_t v;
      // -snan has encoding m=2^63+1, e=65535 (for example)
      v.e = 0xffff;
      v.m = 0x8000000000000001ul;
      *item = v.f;
      (*count)++;
    }
    else if (sscanf(buf, "%La", item) == 1)
    {
      (*count)++;
    }
  }
}

/* check scaled worst-cases from log2l.wc */
static void
check_scaled_worst_cases (void)
{
  long double *items;
  int count;
  readstdin (&items, &count);
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (int i = 0; i < count; i++) {
    ref_init();
    ref_fesetround(rnd);
    fesetround(rnd1[rnd]);
    long double x = items[i];
    check_extended (x);
  }
  free (items);
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

  printf ("   Checking scaled worst cases\n");
  check_scaled_worst_cases ();

  printf ("   Checking random values\n");
#define N 10000000UL /* total number of tests */

  unsigned int seed = getpid ();
  srand (seed);

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

/* Check correctness of bivariate long double function on worst cases.

Copyright (c) 2022-2024 Stéphane Glondu, Paul Zimmermann, Inria.

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

#define _POSIX_C_SOURCE 200809L  /* for getline */

#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#include <mpfr.h>
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif

#include "function_under_test.h"

#include "cm_types.h"

long double cr_function_under_test (long double, long double);
long double ref_function_under_test (long double, long double);
int ref_fesetround (int);
void ref_init (void);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd;

typedef long double ldouble2[2];

/* scanf %La from buf, allowing snan, +snan and -snan */
static int
sscanf_snan (char *buf, long double *x)
{
  if (sscanf(buf, "%La", x) == 1)
    return 1;
  else if (strncmp (buf, "snan", 4) == 0 || strncmp (buf, "+snan", 5) == 0)
  {
    b80u80_t v;
    // +snan has encoding m=2^63+1, e=32767 (for example)
    v.e = 0x7fff;
    v.m = 0x8000000000000001ul;
    *x = v.f;
    return 1;
  }
  else if (strncmp (buf, "-snan", 5) == 0)
  {
    b80u80_t v;
    // -snan has encoding m=2^63+1, e=65535 (for example)
    v.e = 0xffff;
    v.m = 0x8000000000000001ul;
    *x = v.f;
    return 1;
  }
  return 0;
}

static void
readstdin(ldouble2 **result, int *count)
{
  char *buf = NULL;
  size_t buflength = 0;
  ssize_t n;
  int allocated = 512;

  *count = 0;
  if (NULL == (*result = malloc(allocated * sizeof(ldouble2)))) {
    fprintf(stderr, "malloc failed\n");
    exit(1);
  }

  while ((n = getline(&buf, &buflength, stdin)) >= 0) {
    if (n > 0 && buf[0] == '#') continue;
    if (*count >= allocated) {
      int newsize = 2 * allocated;
      ldouble2 *newresult = realloc(*result, newsize * sizeof(ldouble2));
      if (NULL == newresult) {
        fprintf(stderr, "realloc(%d) failed\n", newsize);
        exit(1);
      }
      allocated = newsize;
      *result = newresult;
    }
    ldouble2 *item = *result + *count;
    if (sscanf(buf, "%La,%La", &(*item)[0], &(*item)[1]) == 2)
      (*count)++;
    else if (sscanf_snan (buf, &(*item)[0]) == 1)
    {
      char *tbuf = buf;
      while (*tbuf++ != ',');
      if (sscanf_snan (tbuf, &(*item)[1]) == 1)
        (*count)++;
    }
  }
}

static int
is_nan (long double x)
{
  b80u80_t v = {.f = x};
  return ((v.e&0x7fff) == 0x7fff && (v.m != (1ul << 63)));
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

int tests = 0;

static void
check (long double x, long double y)
{
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp atomic update
#endif
  tests ++;
  ref_init();
  ref_fesetround(rnd);
  mpfr_flags_clear (MPFR_FLAGS_INEXACT);
  long double z1 = ref_function_under_test(x, y);
  mpfr_flags_t inex1 = mpfr_flags_test (MPFR_FLAGS_INEXACT);
  fesetround(rnd1[rnd]);
  feclearexcept (FE_INEXACT);
  long double z2 = cr_function_under_test(x, y);
  fexcept_t inex2;
  fegetexceptflag (&inex2, FE_INEXACT);
  /* Note: the test z1 != z2 would not distinguish +0 and -0. */
  if (is_equal (z1, z2) == 0) {
#ifndef EXCHANGE_X_Y
    printf("FAIL x=%La y=%La ref=%La z=%La\n", x, y, z1, z2);
#else
    printf("FAIL y=%La x=%La ref=%La z=%La\n", x, y, z1, z2);
#endif
    fflush(stdout);
    exit(1);
  }
  if ((inex1 == 0) && (inex2 != 0))
  {
    printf ("Spurious inexact exception for x=%La y=%La (z=%La)\n", x, y, z1);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
  if ((inex1 != 0) && (inex2 == 0))
  {
    printf ("Missing inexact exception for x=%La y=%La (z=%La)\n", x, y, z1);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
}

void
doloop(void)
{
  ldouble2 *items;
  int count;

  readstdin(&items, &count);

#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (int i = 0; i < count; i++) {
    long double x = items[i][0], y = items[i][1];
    check (x, y);
#ifdef WORST_SYMMETRIC_Y
    check (x, -y);
#endif
#ifdef WORST_SYMMETRIC_X
    check (-x, y);
#ifdef WORST_SYMMETRIC_Y
    check (-x, -y);
#endif
#endif
#ifdef WORST_SWAP
    check (y, x);
#ifdef WORST_SYMMETRIC_Y
    check (-y, x);
#endif
#ifdef WORST_SYMMETRIC_X
    check (y, -x);
#ifdef WORST_SYMMETRIC_Y
    check (-y, -x);
#endif
#endif
#endif
  }

  free(items);
  printf("%d tests passed\n", tests);
}

// When x is a NaN, returns 1 if x is an sNaN and 0 if it is a qNaN
static inline int issignaling(long double x) {
  b80u80_t u = {.f = x};

  return !(u.m >> 63);
}

/* check for signaling NaN input */
static void
check_signaling_nan (void)
{
  b80u80_t u;
  u.e = 0x7fffu;
  u.m = 0xc000000000000000ull;
  long double snan = u.f;
  long double y = cr_function_under_test (snan, 2.0L);
  // check that foo(NaN) = NaN
  if (!is_nan (y))
  {
    fprintf (stderr, "Error, foo(sNaN,x) should be NaN, got %La\n", y);
    exit (1);
  }
  // check that the signaling bit disappeared
  if (issignaling (y))
  {
    fprintf (stderr, "Error, foo(sNaN,x) should be qNaN, got %La\n", y);
    exit (1);
  }
  y = cr_function_under_test (2.0L, snan);
  // check that foo(NaN) = NaN
  if (!is_nan (y))
  {
    fprintf (stderr, "Error, foo(x,sNaN) should be NaN, got %La\n", y);
    exit (1);
  }
  // check that the signaling bit disappeared
  if (issignaling (y))
  {
    fprintf (stderr, "Error, foo(x,sNaN) should be qNaN, got %La\n", y);
    exit (1);
  }
  // check also sNaN with the sign bit set
  u.e = 0xffffu;
  u.m = 0xc000000000000000ull;
  snan = u.f;
  y = cr_function_under_test (snan, 2.0L);
  // check that foo(NaN) = NaN
  if (!is_nan (y))
  {
    fprintf (stderr, "Error, foo(sNaN,x) should be NaN, got %La\n", y);
    exit (1);
  }
  // check that the signaling bit disappeared
  if (issignaling (y))
  {
    fprintf (stderr, "Error, foo(sNaN,x) should be qNaN, got %La\n", y);
    exit (1);
  }
  y = cr_function_under_test (2.0L, snan);
  // check that foo(NaN) = NaN
  if (!is_nan (y))
  {
    fprintf (stderr, "Error, foo(x,sNaN) should be NaN, got %La\n", y);
    exit (1);
  }
  // check that the signaling bit disappeared
  if (issignaling (y))
  {
    fprintf (stderr, "Error, foo(x,sNaN) should be qNaN, got %La\n", y);
    exit (1);
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
      else
        {
          fprintf (stderr, "Error, unknown option %s\n", argv[1]);
          exit (1);
        }
    }

  check_signaling_nan ();

  doloop();
}

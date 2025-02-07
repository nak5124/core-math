/* Check correctness of bivariate binary64 function on worst cases.

Copyright (c) 2022-2025 Stéphane Glondu, Paul Zimmermann, Inria.

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
#include <inttypes.h>
#include <fenv.h>
#include <mpfr.h>
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif
#include <errno.h>

#include "function_under_test.h"

double cr_function_under_test (double, double);
double ref_function_under_test (double, double);
int ref_fesetround (int);
void ref_init (void);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd;

typedef double double2[2];

typedef struct {
  double x;
  double y;
#ifdef CORE_MATH_SUPPORT_ERRNO
  int errno_ref;
#endif
} testcase;

typedef union { double f; uint64_t i; } d64u64;

/* scanf %la from buf, allowing snan, +snan and -snan */
static int
sscanf_snan (char *buf, double *x)
{
  if (sscanf(buf, "%la", x) == 1)
    return 1;
  else if (strncmp (buf, "snan", 4) == 0 || strncmp (buf, "+snan", 5) == 0)
  {
    d64u64 u = {.i = 0x7ff4000000000000};
    *x = u.f;
    return 1;
  }
  else if (strncmp (buf, "-snan", 5) == 0)
  {
    d64u64 u = {.i = 0xfff4000000000000};
    *x = u.f;
    return 1;
  }
  return 0;
}

static void
readstdin(testcase **result, int *count)
{
  char *buf = NULL;
  size_t buflength = 0;
  ssize_t n;
  int allocated = 512;

  *count = 0;
  if (NULL == (*result = malloc(allocated * sizeof(testcase)))) {
    fprintf(stderr, "malloc failed\n");
    exit(1);
  }

  while ((n = getline(&buf, &buflength, stdin)) >= 0) {
    if (n > 0 && buf[0] == '#') continue;
    if (*count >= allocated) {
      int newsize = 2 * allocated;
      testcase *newresult = realloc(*result, newsize * sizeof(testcase));
      if (NULL == newresult) {
        fprintf(stderr, "realloc(%d) failed\n", newsize);
        exit(1);
      }
      allocated = newsize;
      *result = newresult;
    }
    testcase *item = *result + *count;
#ifndef CORE_MATH_SUPPORT_ERRNO
    if (sscanf(buf, "%la,%la", &item->x, &item->y) == 2)
      (*count)++;
#else
    char err_str[8];
    int readcnt = sscanf(buf, "%la,%la,%7s", &item->x, &item->y, err_str);
    if (readcnt == 2) {
      item->errno_ref = 0;
      (*count)++;
    }
    else if (readcnt == 3) {
      if (strncmp(err_str, "ERANGE", 7U) == 0) {
        item->errno_ref = ERANGE;
      }
      else if (strncmp(err_str, "EDOM", 5U) == 0) {
        item->errno_ref = EDOM;
      }
      else {
        item->errno_ref = 0;
      }

      (*count)++;
    }
#endif
    else if (sscanf_snan (buf, &item->x) == 1)
    {
      char *tbuf = buf;
      while (*tbuf++ != ',');
      if (sscanf_snan (tbuf, &item->y) == 1)
        (*count)++;
    }
  }
}

static inline uint64_t
asuint64 (double f)
{
  d64u64 u = {.f = f};
  return u.i;
}

/* define our own is_nan function to avoid depending from math.h */
static inline int
is_nan (double x)
{
  uint64_t u = asuint64 (x);
  uint64_t e = u >> 52;
  return (e == 0x7ff || e == 0xfff) && (u << 12) != 0;
}

static inline int
is_equal (double x, double y)
{
  if (is_nan (x))
    return is_nan (y);
  if (is_nan (y))
    return is_nan (x);
  return asuint64 (x) == asuint64 (y);
}

int tests = 0;

// prints snan/qnan when x is NaN
static void
print_binary64 (double x)
{
  d64u64 v = {.f = x};
  if (!is_nan (x)) // not NaN
    printf ("%la", x);
  else
  {
    // if bit 51 is 1, this is a qNaN, otherwise a sNaN
    if ((v.i >> 51) & 1)
      printf ("qnan");
    else
      printf ("snan");
  }
}

int underflow_before; // non-zero if processor raises underflow before rounding

// return non-zero if the processor raises underflow before rounding
// (e.g., aarch64)
static void
check_underflow_before (void)
{
  fexcept_t flag;
  fegetexceptflag (&flag, FE_ALL_EXCEPT); // save flags
  fesetround (FE_TONEAREST);
  feclearexcept (FE_UNDERFLOW);
  float x = 0x1p-126f;
  float y = __builtin_fmaf (-x, x, x);
  if (x == y) // this is needed otherwise the compiler says y is unused
    underflow_before = fetestexcept (FE_UNDERFLOW);
  fesetexceptflag (&flag, FE_ALL_EXCEPT); //restore flags
}

/* In case of underflow before rounding and |z| = 2^-1022, raises the MPFR
   underflow exception if |f(x,y)| < 2^-1022. */
static void
fix_spurious_underflow (double x, double y, double z)
{
  if (!underflow_before || __builtin_fabs (z) != 0x1p-1022)
    return;
  // the processor raises underflow before rounding, and |z| = 2^-1022
  mpfr_t t, u;
  mpfr_init2 (t, 53);
  mpfr_init2 (u, 53);
  mpfr_set_d (t, x, MPFR_RNDN); // exact
  mpfr_set_d (u, y, MPFR_RNDN); // exact
  mpfr_function_under_test (t, t, u, MPFR_RNDZ);
  mpfr_abs (t, t, MPFR_RNDN); // exact
  if (mpfr_cmp_d (t, 0x1p-1022) < 0) // |f(x,y)| < 2^-1022
    mpfr_set_underflow ();
  mpfr_clear (t);
  mpfr_clear (u);
}

static void
check (testcase ts)
{
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp atomic update
#endif
  tests ++;
  ref_init();
  ref_fesetround(rnd);
  mpfr_flags_clear (MPFR_FLAGS_INEXACT);
  double z1 = ref_function_under_test(ts.x, ts.y);
#ifdef CORE_MATH_CHECK_INEXACT
  mpfr_flags_t inex1 = mpfr_flags_test (MPFR_FLAGS_INEXACT);
#endif
  fesetround(rnd1[rnd]);
  feclearexcept (FE_INEXACT);
#ifdef CORE_MATH_SUPPORT_ERRNO
  errno = 0;
#endif
  double z2 = cr_function_under_test(ts.x, ts.y);
#ifdef CORE_MATH_SUPPORT_ERRNO
  int cr_errno = errno;
#endif
  /* Note: the test z1 != z2 would not distinguish +0 and -0. */
  if (is_equal (z1, z2) == 0) {
#ifndef EXCHANGE_X_Y
    printf("FAIL x=");
    print_binary64 (ts.x);
    printf (" y=");
    print_binary64 (ts.y);
    printf (" ref=");
    print_binary64 (z1);
    printf (" z=");
    print_binary64 (z2);
    printf ("\n");
#else
    printf("FAIL y=");
    print_binary64 (ts.x);
    printf (" x=");
    print_binary64 (ts.y);
    printf (" ref=");
    print_binary64 (z1);
    printf (" z=");
    print_binary64 (z2);
    printf ("\n");
#endif
    fflush(stdout);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }

  fix_spurious_underflow (ts.x, ts.y, z1);

  // Check for spurious/missing underflow exception
  if (fetestexcept (FE_UNDERFLOW) && !mpfr_flags_test (MPFR_FLAGS_UNDERFLOW)
      && !is_nan (z1))
  {
    printf ("Spurious underflow exception for x,y=%la,%la (z=%la)\n",
            ts.x, ts.y, z1);
    fflush (stdout);
#ifdef DO_NOT_ABORT
    return 1;
#else
    exit(1);
#endif
  }
  if (!fetestexcept (FE_UNDERFLOW) && mpfr_flags_test (MPFR_FLAGS_UNDERFLOW))
  {
    printf ("Missing underflow exception for x,y=%la,%la (z=%la)\n",
            ts.x, ts.y, z1);
    fflush (stdout);
#ifdef DO_NOT_ABORT
    return 1;
#else
    exit(1);
#endif
  }

  /* Check for spurious/missing overflow exception */
  if (fetestexcept (FE_OVERFLOW) && !mpfr_flags_test (MPFR_FLAGS_OVERFLOW)
      && !is_nan (z1))
  {
    printf ("Spurious overflow exception for x,y=%la,%la (z=%la)\n",
            ts.x, ts.y, z1);
    fflush (stdout);
#ifdef DO_NOT_ABORT
    return 1;
#else
    exit(1);
#endif
  }
  if (!fetestexcept (FE_OVERFLOW) && mpfr_flags_test (MPFR_FLAGS_OVERFLOW))
  {
    printf ("Missing overflow exception for x,y=%la,%la (z=%la)\n",
            ts.x, ts.y, z1);
    fflush (stdout);
#ifdef DO_NOT_ABORT
    return 1;
#else
    exit(1);
#endif
  }

#ifdef CORE_MATH_CHECK_INEXACT
  int inex2 = fetestexcept (FE_INEXACT);
  if ((inex1 == 0) && (inex2 != 0))
  {
    printf ("Spurious inexact exception for x=%la y=%la (z=%la)\n", ts.x, ts.y, z1);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
  if ((inex1 != 0) && (inex2 == 0))
  {
    printf ("Missing inexact exception for x=%la y=%la (z=%la)\n", ts.x, ts.y, z1);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
#endif

#ifdef CORE_MATH_SUPPORT_ERRNO
  // most tests don't check for errno setting, so it's not yet possible to check when errno was set incorrectly (case when errno_ref = 0 & cr_errno != 0)
  if (ts.errno_ref != 0 && cr_errno != ts.errno_ref) {
    printf("%s error not set for x,y=%la,%la (z=%la)\n", ts.errno_ref == ERANGE ? "Range" : "Domain", ts.x, ts.y, z1);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
#endif
}

void
doloop(void)
{
  testcase *items;
  int count;

  readstdin(&items, &count);

#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (int i = 0; i < count; i++) {
    testcase ts = items[i];
    double x = ts.x, y = ts.y;

    /* the following is just to avoid a compiler warning that x and y are
       unused when neither WORST_SYMMETRIC_* no WORST_SWAP is defined */
    ts.x = x;
    ts.y = y;
    check (ts);
#ifdef WORST_SYMMETRIC_Y
    ts.x = x;
    ts.y = -y;
    check (ts);
#endif
#ifdef WORST_SYMMETRIC_X
    ts.x = -x;
    ts.y = y;
    check (ts);
#ifdef WORST_SYMMETRIC_Y
    ts.x = -x;
    ts.y = -y;
    check (ts);
#endif
#endif
#ifdef WORST_SWAP
    ts.x = y;
    ts.y = x;
    check (ts);
#ifdef WORST_SYMMETRIC_Y
    ts.x = -y;
    ts.y = x;
    check (ts);
#endif
#ifdef WORST_SYMMETRIC_X
    ts.x = y;
    ts.y = -x;
    check (ts);
#ifdef WORST_SYMMETRIC_Y
    ts.x = -y;
    ts.y = -x;
    check (ts);
#endif
#endif
#endif
  }

  free(items);
  printf("%d tests passed\n", tests);
}

static inline double
asfloat64 (uint64_t i)
{
  d64u64 u = {.i = i};
  return u.f;
}

// When x is a NaN, returns 1 if x is an sNaN and 0 if it is a qNaN
static inline int issignaling(double x) {
  d64u64 _x = {.f = x};

  return !(_x.i & (1ull << 51));
}

/* check for signaling NaN input */
static void
check_signaling_nan (void)
{
  double snan = asfloat64 (0x7ff0000000000001ull);
  double y = cr_function_under_test (snan, 1.0);
  // check that foo(NaN) = NaN
  if (!is_nan (y))
  {
    fprintf (stderr, "Error, foo(sNaN,1.0) should be NaN, got %la=%"PRIx64"\n",
             y, asuint64 (y));
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }
  // check that the signaling bit disappeared
  if (issignaling (y))
  {
    fprintf (stderr, "Error, foo(sNaN,1.0) should be qNaN, got sNaN=%"PRIx64"\n",
             asuint64 (y));
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }
  y = cr_function_under_test (1.0, snan);
  // check that foo(NaN) = NaN
  if (!is_nan (y))
  {
    fprintf (stderr, "Error, foo(1.0,sNaN) should be NaN, got %la=%"PRIx64"\n",
             y, asuint64 (y));
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }
  // check that the signaling bit disappeared
  if (issignaling (y))
  {
    fprintf (stderr, "Error, foo(1.0,sNaN) should be qNaN, got sNaN=%"PRIx64"\n",
             asuint64 (y));
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }
  // check also sNaN with the sign bit set
  snan = asfloat64 (0xfff0000000000001ull);
  y = cr_function_under_test (snan, 1.0);
  // check that foo(NaN) = NaN
  if (!is_nan (y))
  {
    fprintf (stderr, "Error, foo(sNaN,1.0) should be NaN, got %la=%"PRIx64"\n",
             y, asuint64 (y));
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }
  // check that the signaling bit disappeared
  if (issignaling (y))
  {
    fprintf (stderr, "Error, foo(sNaN,1.0) should be qNaN, got sNaN=%"PRIx64"\n",
             asuint64 (y));
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }
  y = cr_function_under_test (1.0, snan);
  // check that foo(NaN) = NaN
  if (!is_nan (y))
  {
    fprintf (stderr, "Error, foo(1.0,sNaN) should be NaN, got %la=%"PRIx64"\n",
             y, asuint64 (y));
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }
  // check that the signaling bit disappeared
  if (issignaling (y))
  {
    fprintf (stderr, "Error, foo(1.0,sNaN) should be qNaN, got sNaN=%"PRIx64"\n",
             asuint64 (y));
#ifndef DO_NOT_ABORT
    exit (1);
#endif
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

  check_underflow_before ();

  check_signaling_nan ();

  doloop();
}

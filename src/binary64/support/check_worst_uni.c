/* Check correctness of univariate binary64 function on worst cases.

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
#ifdef CORE_MATH_SUPPORT_ERRNO
#include <errno.h>
#endif

#include "cm_types.h"

#include "function_under_test.h"

double cr_function_under_test (double);
double ref_function_under_test (double);
int ref_fesetround (int);
void ref_init (void);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;

typedef struct {
  double x;
#ifdef CORE_MATH_SUPPORT_ERRNO
  int errno_ref;
#endif
} testcase;

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

    // special code for snan, since glibc does not read them
    if (strncmp (buf, "snan", 4) == 0 || strncmp (buf, "+snan", 5) == 0)
    {
      /* According to IEEE 754-2019, qNaN's have 1 as upper bit of their
         52-bit significand, and sNaN's have 0 */
      b64u64_u u = {.u = 0x7ff4000000000000};
      item->x = u.f;
      (*count)++;
    }
    else if (strncmp (buf, "-snan", 5) == 0)
    {
      b64u64_u u = {.u = 0xfff4000000000000};
      item->x = u.f;
      (*count)++;
    }
#ifndef CORE_MATH_SUPPORT_ERRNO
    else if (sscanf(buf, "%la", &(item->x)) == 1)
      (*count)++;
#else
    else {
      char err_str[8];
      int readcnt = sscanf(buf, "%la,%7s", &(item->x), err_str);
      if (readcnt == 1) {
        item->errno_ref = 0;
        (*count)++;
      }
      else if (readcnt == 2) {
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
    }
#endif
  }
}

static inline uint64_t
asuint64 (double f)
{
  b64u64_u u = {.f = f};
  return u.u;
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

// return 1 if failure, 0 otherwise
static int
check (testcase ts)
{
  ref_init();
  ref_fesetround(rnd);
  mpfr_flags_clear (MPFR_FLAGS_INEXACT);
  double z1 = ref_function_under_test(ts.x);
  mpfr_flags_t inex1 = mpfr_flags_test (MPFR_FLAGS_INEXACT);
  fesetround(rnd1[rnd]);
  feclearexcept (FE_INEXACT);
#ifdef CORE_MATH_SUPPORT_ERRNO
  errno = 0;
#endif
  double z2 = cr_function_under_test(ts.x);
#ifdef CORE_MATH_SUPPORT_ERRNO
  int cr_errno = errno;
#endif
  fexcept_t inex2;
  fegetexceptflag (&inex2, FE_INEXACT);
  /* Note: the test z1 != z2 would not distinguish +0 and -0. */
  if (is_equal (z1, z2) == 0) {
    printf("FAIL x=%la ref=%la z=%la\n", ts.x, z1, z2);
    fflush(stdout);
#ifdef DO_NOT_ABORT
    return 1;
#else
    exit(1);
#endif
  }
#ifdef CORE_MATH_CHECK_INEXACT
  if ((inex1 == 0) && (inex2 != 0))
  {
    printf ("Spurious inexact exception for x=%la (y=%la)\n", ts.x, z1);
    fflush (stdout);
#ifdef DO_NOT_ABORT
    return 1;
#else
    exit(1);
#endif
  }
  if ((inex1 != 0) && (inex2 == 0))
  {
    printf ("Missing inexact exception for x=%la (y=%la)\n", ts.x, z1);
    fflush (stdout);
#ifdef DO_NOT_ABORT
    return 1;
#else
    exit(1);
#endif
  }
#endif

#ifdef CORE_MATH_SUPPORT_ERRNO
  // most tests don't check for errno setting, so it's not yet possible to check when errno was set incorrectly (case when errno_ref = 0 & cr_errno != 0)
  if (ts.errno_ref != 0 && cr_errno != ts.errno_ref) {
    printf("%s error not set for x=%la (y=%la)\n", ts.errno_ref == ERANGE ? "Range" : "Domain", ts.x, z1);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
#endif
  return 0;
}

void
doloop(void)
{
  testcase *items;
  int count, tests = 0, failures = 0;

  readstdin(&items, &count);

#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for reduction(+: failures,tests)
#endif
  for (int i = 0; i < count; i++) {
    testcase ts = items[i];
    tests ++;
    if (check (ts))
      failures ++;
#ifdef WORST_SYMMETRIC
    tests ++;
    ts.x = -ts.x;
    if (check (ts))
      failures ++;
#endif
  }

  free(items);
  printf("%d tests passed, %d failure(s)\n", tests, failures);
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

  doloop();
}

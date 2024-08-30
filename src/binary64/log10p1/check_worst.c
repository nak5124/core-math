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
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif

#include "function_under_test.h"

double cr_function_under_test (double);
double ref_function_under_test (double);
int ref_fesetround (int);
void ref_init (void);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;

typedef union { double f; uint64_t i; } d64u64;

static void
readstdin(double **result, int *count)
{
  char *buf = NULL;
  size_t buflength = 0;
  ssize_t n;
  int allocated = 512;

  *count = 0;
  if (NULL == (*result = malloc(allocated * sizeof(double)))) {
    fprintf(stderr, "malloc failed\n");
    exit(1);
  }

  while ((n = getline(&buf, &buflength, stdin)) >= 0) {
    if (n > 0 && buf[0] == '#') continue;
    if (*count >= allocated) {
      int newsize = 2 * allocated;
      double *newresult = realloc(*result, newsize * sizeof(double));
      if (NULL == newresult) {
        fprintf(stderr, "realloc(%d) failed\n", newsize);
        exit(1);
      }
      allocated = newsize;
      *result = newresult;
    }
    double *item = *result + *count;
    // special code for snan, since glibc does not read them
    if (strncmp (buf, "snan", 4) == 0 || strncmp (buf, "+snan", 5) == 0)
    {
      /* According to IEEE 754-2019, qNaN's have 1 as upper bit of their
         52-bit significand, and sNaN's have 0 */
      d64u64 u = {.i = 0x7ff4000000000000};
      *item = u.f;
      (*count)++;
    }
    else if (strncmp (buf, "-snan", 5) == 0)
    {
      d64u64 u = {.i = 0xfff4000000000000};
      *item = u.f;
      (*count)++;
    }
    else if (sscanf(buf, "%la", item) == 1)
      (*count)++;
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
  int e = u >> 52;
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

int tests = 0, failures = 0;

static void
check (double x)
{
  ref_init();
  ref_fesetround(rnd);
  fesetround(rnd1[rnd]);
  double z1 = ref_function_under_test(x);
  double z2 = cr_function_under_test(x);
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp atomic update
#endif
  tests ++;
  /* Note: the test z1 != z2 would not distinguish +0 and -0. */
  if (is_equal (z1, z2) == 0) {
    printf("FAIL x=%la ref=%la z=%la\n", x, z1, z2);
    fflush(stdout);
#ifdef DO_NOT_ABORT
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp atomic update
#endif
    failures ++;
#else
    exit(1);
#endif
  }
}

void
doloop(void)
{
  double *items;
  int count;

  readstdin(&items, &count);

#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (int i = 0; i < count; i++) {
    double x = items[i];
    check (x);
#ifdef WORST_SYMMETRIC
    check (-x);
#endif
    // if x in [2^-1020, 2^-1019), check all 2^-1019 <= 2^k*x < 2^-109
    if (0x1p-1020 <= x && x < 0x1p-1019)
    {
      while (x < 0x1p-110)
      {
        x = 2.0 * x;
        check (x);
#ifdef WORST_SYMMETRIC
        check (-x);
#endif
      }
    }
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

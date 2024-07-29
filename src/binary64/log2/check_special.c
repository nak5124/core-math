/* Generate special cases for log2 testing.

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
#include <assert.h>

int ref_init (void);
int ref_fesetround (int);

double cr_log2 (double);
double ref_log2 (double);

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

static inline double
asfloat64 (uint64_t i)
{
  union
  {
    uint64_t i;
    double f;
  } u = {i};
  return u.f;
}

static void
check (double x)
{
  double y1 = ref_log2 (x);
  fesetround (rnd1[rnd]);
  double y2 = cr_log2 (x);
  if (isnan (y1) && isnan (y2))
    return;
  if (asuint64 (y1) != asuint64 (y2))
  {
    printf ("FAIL x=%la ref=%la z=%la\n", x, y1, y2);
    fflush (stdout);
    exit (1);
  }
}

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

/* check scaled worst-cases from log2.wc */
static void
check_scaled_worst_cases (void)
{
  double *items;
  int count, tests, failures;
  readstdin (&items, &count);
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for reduction(+: failures,tests)
#endif
  for (int i = 0; i < count; i++) {
    ref_init();
    ref_fesetround(rnd);
    fesetround(rnd1[rnd]);
    double x1 = items[i];
    if (!isnan (x1) && 2 * x1 != x1) // is is not NaN nor +/Inf nor +/-0
      {
        int e;
        double x0 = frexp (x1, &e);
        for (e = -1074; e <= 1024; e++)
          {
            double x = ldexp (x0, e);
            double z1 = ref_log2 (x);
            double z2 = cr_log2 (x);
            tests ++;
            /* Note: the test z1 != z2 would not distinguish +0 and -0. */
            if (z1 != z2) {
              printf("FAIL x1=%la x=%la ref=%la z=%la\n", x1, x, z1, z2);
              fflush(stdout);
#ifdef DO_NOT_ABORT
              failures ++;
#else
              exit(1);
#endif
            }
        }
    }
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
  ref_init ();
  ref_fesetround (rnd);

  printf ("Checking scaled worst cases...\n");
  check_scaled_worst_cases ();

#define K 1000000000UL /* total number of tests */
#define BUF_SIZE 1000

  long seed = getpid ();
  srand48 (seed);
  
  double buf[BUF_SIZE];
  uint64_t N = K / BUF_SIZE;
  printf ("Checking random numbers...\n");
  for (uint64_t n = 0; n < N; n++)
  {
    /* warning: lrand48 is not thread-safe, thus we put it outside
       the parallel loop */
    for (int i = 0; i < BUF_SIZE; i++)
    {
      uint64_t j = ((uint64_t) lrand48 () << 62)
        | ((uint64_t) lrand48 () << 31) | (uint64_t) lrand48 ();
      double x = asfloat64 (j);
      buf[i] = (x >= 0) ? x : -x;
    }
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
    for (int i = 0; i < BUF_SIZE; i++)
      check (buf[i]);
  }

  return 0;
}

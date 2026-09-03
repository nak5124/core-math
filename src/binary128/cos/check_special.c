/* Generate special cases for cosq testing.

Copyright (c) 2026 Chen-Pang He <jdh863@gmail.com>

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
#include <inttypes.h>
#include <string.h>
#include <fenv.h>
#include <math.h>
#include <unistd.h>
#include <quadmath.h>
#ifdef CORE_MATH_SUPPORT_ERRNO
#include <errno.h>
#endif
#ifndef MPFR_WANT_FLOAT128
#define MPFR_WANT_FLOAT128
#endif
#include <mpfr.h>
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif

typedef unsigned __int128 u128;

typedef union {
  u128 a;
  __float128 f;
} b128u128_u;

int ref_fesetround(int);
void ref_init(void);

__float128 cr_cosq(__float128);
__float128 ref_cosq(__float128);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

#define MAX_THREADS 192

static unsigned int Seed[MAX_THREADS];

#ifndef CORE_MATH_TESTS
/* total number of tests (mpfr_cos is slow for huge exponents) */
#define CORE_MATH_TESTS (100UL*1000*1000)
#endif

/* define our own is_nan function to avoid depending from math.h */
static inline int is_nan(__float128 x){
  b128u128_u u = {.f = x};
  return (u.a<<1)>((u128)0x7fff<<113);
}

/* when x is a NaN, returns 1 if x is an sNaN and 0 if it is a qNaN */
static inline int is_signaling(__float128 x){
  b128u128_u u = {.f = x};
  return !(u.a & ((u128)1<<111));
}

static inline int is_equal(__float128 x, __float128 y){
  if(is_nan(x)) return is_nan(y);
  if(is_nan(y)) return is_nan(x);
  b128u128_u ux = {.f = x}, uy = {.f = y};
  return ux.a == uy.a;
}

/* We have to write some special code to print errors since:
   (1) the printf hook for quadmath is not installed everywhere,
       thus printf ("%Qa", ...) does not work everywhere
   (2) quadmath_snprintf can only print a single quad,
       but quadmath_snprintf (..., "x=%Qa", x) does not work
*/
static void
error2 (__float128 x, __float128 y, __float128 z)
{
  char buf[256];
  int n;
  n = snprintf (buf, sizeof buf, "FAIL x=");
  n += quadmath_snprintf (buf + n, sizeof buf - n, "%Qa", x);
  n += snprintf (buf + n, sizeof buf - n, " ref=");
  n += quadmath_snprintf (buf + n, sizeof buf - n, "%Qa", y);
  n += snprintf (buf + n, sizeof buf - n, " z=");
  n += quadmath_snprintf (buf + n, sizeof buf - n, "%Qa", z);
  printf ("%s\n", buf);
}

/* print "FAIL <err> for x=<x> (y=<y>)", with the bit patterns of x and y
   since %Qa does not show the payload of a NaN */
static void
error (const char *err, __float128 x, __float128 y)
{
  b128u128_u ux = {.f = x}, uy = {.f = y};
  char buf[256];
  int n;
  n = snprintf (buf, sizeof buf, "FAIL %s for x=", err);
  n += quadmath_snprintf (buf + n, sizeof buf - n, "%Qa", x);
  n += snprintf (buf + n, sizeof buf - n, " [0x%016" PRIx64 "%016" PRIx64
                 "] (y=", (uint64_t) (ux.a >> 64), (uint64_t) ux.a);
  n += quadmath_snprintf (buf + n, sizeof buf - n, "%Qa", y);
  n += snprintf (buf + n, sizeof buf - n, " [0x%016" PRIx64 "%016" PRIx64
                 "])", (uint64_t) (uy.a >> 64), (uint64_t) uy.a);
  printf ("%s\n", buf);
  fflush (stdout);
}

static void check(__float128 x){
  ref_init();
  ref_fesetround(rnd);
  __float128 y1 = ref_cosq(x), y2 = cr_cosq(x);
  if(!is_equal(y1, y2)) {
    error2 (x, y1, y2);
    fflush(stdout);
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }
}

/* check that cr_cosq(x) is a qNaN, that the invalid exception is raised
   iff expect_invalid, and (with CORE_MATH_SUPPORT_ERRNO) that errno is set
   to EDOM iff expect_edom, and left unchanged otherwise */
static void
check_nan (__float128 x, int expect_invalid, int expect_edom)
{
  (void) expect_edom; /* only used with CORE_MATH_SUPPORT_ERRNO */
  feclearexcept (FE_ALL_EXCEPT);
#ifdef CORE_MATH_SUPPORT_ERRNO
  errno = 0;
#endif
  __float128 y = cr_cosq (x);
  if (!is_nan (y))
    error ("cosq(x) should be NaN", x, y);
  else if (is_signaling (y))
    error ("cosq(x) should be a qNaN", x, y);
  else if (expect_invalid && !fetestexcept (FE_INVALID))
    error ("Missing invalid exception", x, y);
  else if (!expect_invalid && fetestexcept (FE_INVALID))
    error ("Spurious invalid exception", x, y);
#ifdef CORE_MATH_SUPPORT_ERRNO
  else if (expect_edom && errno != EDOM)
    error ("Missing errno=EDOM", x, y);
  else if (!expect_edom && errno != 0)
    error ("Spurious errno", x, y);
#endif
  else
    return;
#ifndef DO_NOT_ABORT
  exit (1);
#endif
}

/* check that cr_cosq(x) equals expected bit-for-bit, raises no exception
   (in particular no inexact exception), and (with CORE_MATH_SUPPORT_ERRNO)
   leaves errno unchanged */
static void
check_exact (__float128 x, __float128 expected)
{
  feclearexcept (FE_ALL_EXCEPT);
#ifdef CORE_MATH_SUPPORT_ERRNO
  errno = 0;
#endif
  __float128 y = cr_cosq (x);
  b128u128_u uy = {.f = y}, ue = {.f = expected};
  if (uy.a != ue.a)
    error2 (x, expected, y);
  else if (fetestexcept (FE_INEXACT))
    error ("Spurious inexact exception", x, y);
  else if (fetestexcept (FE_ALL_EXCEPT))
    error ("Spurious exception", x, y);
#ifdef CORE_MATH_SUPPORT_ERRNO
  else if (errno != 0)
    error ("Spurious errno", x, y);
#endif
  else
    return;
#ifndef DO_NOT_ABORT
  exit (1);
#endif
}

static void
check_invalid (void)
{
  b128u128_u plus_inf = {.a = (u128) 0x7fff << 112};
  b128u128_u minus_inf = {.a = (u128) 0xffff << 112};
  b128u128_u plus_qnan = {.a = ((u128) 0x7fff << 112) | ((u128) 1 << 111)};
  b128u128_u minus_qnan = {.a = ((u128) 0xffff << 112) | ((u128) 1 << 111)};
  b128u128_u plus_snan = {.a = ((u128) 0x7fff << 112) | 1};
  b128u128_u minus_snan = {.a = ((u128) 0xffff << 112) | 1};
  b128u128_u plus_zero = {.a = 0};
  b128u128_u minus_zero = {.a = (u128) 1 << 127};
  b128u128_u one = {.a = (u128) 0x3fff << 112};

  /* cosq(+/-Inf) = NaN with invalid exception and errno=EDOM */
  check_nan (plus_inf.f, 1, 1);
  check_nan (minus_inf.f, 1, 1);
  /* cosq(qNaN) = qNaN with no exception (and errno unchanged) */
  check_nan (plus_qnan.f, 0, 0);
  check_nan (minus_qnan.f, 0, 0);
  /* cosq(sNaN) = qNaN with invalid exception (and errno unchanged) */
  check_nan (plus_snan.f, 1, 0);
  check_nan (minus_snan.f, 1, 0);
  /* cosq(+0) = cosq(-0) = +1 exactly, with no exception */
  check_exact (plus_zero.f, one.f);
  check_exact (minus_zero.f, one.f);
}

/* random 64-bit value from three rand_r() calls (31 random bits each) */
static uint64_t rand64(int tid){
  uint64_t r = rand_r(Seed + tid);
  r = (r << 31) | rand_r(Seed + tid);
  r = (r << 2) | (rand_r(Seed + tid) & 3);
  return r;
}

static u128 rand128(int tid){
  u128 r = rand64(tid);
  return (r << 64) | rand64(tid);
}

/* random bit pattern (with random sign): this covers all exponents,
   including subnormals, Inf and NaN */
static __float128 get_random(int tid){
  b128u128_u v = {.a = rand128(tid)};
  return v.f;
}

/* random significand and sign, with unbiased exponent uniform in [-60,20]
   (the small-argument and argument-reduction bands) */
static __float128 get_random_banded(int tid){
  b128u128_u v = {.a = rand128(tid)};
  int e = -60 + (int) (rand_r(Seed + tid) % 81);
  v.a &= ~((u128) 0x7fff << 112);
  v.a |= (u128) (e + 16383) << 112;
  return v.f;
}

/* return k*pi/2 rounded to binary128, where k is a random integer with a
   random bit-length, perturbed by -3..+3 ulps, with a random sign: this
   exercises the cancellation in the argument reduction near multiples of
   pi/2. Half of the draws take k < 2^112, where the rounding error of x is
   much smaller than pi/2 (thus x is really within a few ulps of a multiple
   of pi/2); the other half spread the bit-length over [1,16388], which
   covers all binades up to the overflow threshold. */
static __float128 get_random_near_pio2(int tid){
  unsigned long nbits = (rand_r(Seed + tid) & 1)
    ? 1 + rand_r(Seed + tid) % 112
    : 1 + rand_r(Seed + tid) % (16383 + 5);
  unsigned long chunks = (nbits + 63) / 64;
  unsigned long top = nbits - 64 * (chunks - 1); /* bits in the top chunk */
  /* k*pi/2 can exceed the binary128 range: widen the MPFR exponent range
     for this computation (check() restores it through ref_init()) */
  mpfr_exp_t emax = mpfr_get_emax ();
  mpfr_set_emax (mpfr_get_emax_max ());
  mpfr_t k, p;
  mpfr_init2 (k, nbits);
  mpfr_set_ui (k, 0, MPFR_RNDN);
  for (unsigned long i = 0; i < chunks; i++) {
    uint64_t c = rand64 (tid);
    if (i == 0) {
      if (top < 64)
        c &= ((uint64_t) 1 << top) - 1;
      c |= (uint64_t) 1 << (top - 1); /* force the bit-length of k */
    }
    /* k = 2^64*k + c, exact since k has nbits bits of precision */
    mpfr_mul_2ui (k, k, 32, MPFR_RNDN);
    mpfr_add_ui (k, k, (unsigned long) (c >> 32), MPFR_RNDN);
    mpfr_mul_2ui (k, k, 32, MPFR_RNDN);
    mpfr_add_ui (k, k, (unsigned long) (c & 0xffffffff), MPFR_RNDN);
  }
  mpfr_init2 (p, nbits + 130);
  mpfr_const_pi (p, MPFR_RNDN);
  mpfr_mul (p, p, k, MPFR_RNDN);
  mpfr_div_2ui (p, p, 1, MPFR_RNDN);
  b128u128_u v;
  if (mpfr_get_exp (p) > 16384) /* k*pi/2 >= 2^16384 overflows */
    v.a = (u128) 0x7fff << 112;
  else
    v.f = mpfr_get_float128 (p, MPFR_RNDN);
  mpfr_clear (k);
  mpfr_clear (p);
  mpfr_set_emax (emax);
  /* perturb by -3..+3 ulps (p > 0 thus v.f > 0, and the bit pattern of a
     positive binary128 is monotonic in its value) */
  v.a += (u128) (__int128) ((int) (rand_r (Seed + tid) % 7) - 3);
  /* random sign */
  if (rand_r (Seed + tid) & 1)
    v.a ^= (u128) 1 << 127;
  return v.f;
}

static void
check_random (__float128 (*get)(int))
{
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for(uint64_t n = 0; n < CORE_MATH_TESTS / 3; n++){
    ref_init();
    ref_fesetround(rnd);
    int tid;
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    __float128 x = get(tid);
    check(x);
  }
}

int main(int argc, char *argv[]){
  while(argc >= 2){
    if(strcmp(argv[1], "--rndn") == 0){
      rnd = 0;
      argc --;
      argv ++;
    } else if(strcmp(argv[1], "--rndz") == 0){
      rnd = 1;
      argc --;
      argv ++;
    } else if(strcmp(argv[1], "--rndu") == 0){
      rnd = 2;
      argc --;
      argv ++;
    } else if(strcmp(argv[1], "--rndd") == 0){
      rnd = 3;
      argc --;
      argv ++;
    } else if(strcmp(argv[1], "--verbose") == 0){
      verbose = 1;
      argc --;
      argv ++;
    } else {
      fprintf(stderr, "Error, unknown option %s\n", argv[1]);
      exit(1);
    }
  }

  ref_init();
  ref_fesetround(rnd);
  fesetround(rnd1[rnd]);

  unsigned int seed = getpid();
  for(int i = 0; i < MAX_THREADS; i++)
    Seed[i] = seed + i;

  printf("Checking special values\n");
  check_invalid ();

  printf("Checking random values in the full range\n");
  check_random (get_random);

  printf("Checking random values with exponent in [-60,20]\n");
  check_random (get_random_banded);

  printf("Checking random values near multiples of pi/2\n");
  check_random (get_random_near_pio2);

  return 0;
}

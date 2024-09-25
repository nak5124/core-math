#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#include <mpfr.h>
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif
#include <unistd.h>
#include <math.h>
#include <assert.h>

void doloop (int, int);
extern long double cr_powl (long double, long double);
extern int ref_fesetround (int);
extern void ref_init (void);
extern mpfr_rnd_t rnd2[];
extern long double ref_powl (long double, long double);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd;
int verbose = 0;

typedef union {long double f; struct {uint64_t m; uint16_t e;};} b80u80_t;

static long double
get_random ()
{
	b80u80_t v;
  v.m = rand ();
  v.m |= (uint64_t) rand () << 32;
  v.e = rand () & 0xffff;
  // If v is not a denormal, m should have its msb set,
	// otherwise it should be cleared
  uint64_t t = (v.e&0x7fff) != 0;
  v.m = (t << 63) | (v.m & ~(1ul<<63));
  return v.f;
}

/* define our own is_nan function to avoid depending from math.h */
static inline int
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
  return v.e == w.e && v.m == w.m;
}

// return 1 in case of failure, 0 otherwise
static int
check_aux (long double x, long double y)
{
  long double z, t;
  mpfr_t X, Y, Z;
  ref_init();
  ref_fesetround(rnd);
  mpfr_init2 (X, 64);
  mpfr_init2 (Y, 64);
  mpfr_init2 (Z, 64);
  mpfr_set_ld (X, x, MPFR_RNDN);
  mpfr_set_ld (Y, y, MPFR_RNDN);
  mpfr_flags_clear (MPFR_FLAGS_INEXACT);
  t = ref_powl (x, y);
  mpfr_flags_t inex1 = mpfr_flags_test (MPFR_FLAGS_INEXACT);
  fesetround(rnd1[rnd]);
  feclearexcept (FE_INEXACT);
  z = cr_powl (x, y);
  fexcept_t inex2;
  fegetexceptflag (&inex2, FE_INEXACT);
  mpfr_clear (X);
  mpfr_clear (Y);
  mpfr_clear (Z);
  if (!is_equal (z, t))
  {
    printf("FAIL x,y=%La,%La ref=%La z=%La\n", x,y,t,z);
#ifdef DO_NOT_ABORT
    return 1;
#else
    exit(1);
#endif
  }
  if ((inex1 == 0) && (inex2 != 0))
  {
    printf ("Spurious inexact exception for x,y=%La,%La (z=%La)\n", x, y, z);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
  if ((inex1 != 0) && (inex2 == 0))
  {
    printf ("Missing inexact exception for x,y=%La,%La (z=%La)\n", x, y, z);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
  return 0;
}

// return the number of failures
static int
check (long double x, long double y)
{
  int ret = check_aux (x, y);
  // if y is an integer, also check with -x
  if (y == (long double) (int64_t) y)
    ret += check_aux (-x, y);
  return ret;
}

// check x=2^n and y, return 1 iff x is in the long double range
static void
check_pow2_aux (int n, long double y)
{
  if (n < -16445 || 16384 <= n)
    return;
  long double x = ldexpl (1.0L, n);
  check (x, y);
}

// check exact values for x=2^n and y=m/2^k with m odd, k >= 6
static void
check_pow2 (void)
{
  // since n should be a multiple of 2^k and n <= 16445, we have k <= 14
  for (int k = 6; k <= 14; k++)
  {
    int K = 1 << k;
    // positive n
    for (int n = K; n <= 16445; n += K)
    {
      int e = n / K;
      for (int m = 1; m * e <= 16445; m += 2)
      {
        check_pow2_aux (n, ldexpl ((long double) m, -k));
        check_pow2_aux (-n, -ldexpl ((long double) m, -k));
      }
    }
    // negative n
    for (int n = -K; n >= -16445; n -= K)
    {
      int e = n / K;
      for (int m = 1; m * (-e) <= 16445; m += 2)
      {
        check_pow2_aux (n, ldexpl ((long double) m, -k));
        check_pow2_aux (-n, -ldexpl ((long double) m, -k));
      }
    }
  }
}

// perform N random tests near underflow threshold
static void
check_near_underflow (int N)
{
  long double threshold1 = -16446.0L; // half smallest subnormal
  long double threshold2 = -16445.0L; // smallest subnormal
  long double threshold3 = -16382.0L; // smallest normal
  for (int n = 0; n < N / 3; n++)
  {
    long double x = get_random ();
    x = fabsl (x);
    long double y = threshold1 / log2l (x);
    check (x, y);
    y = threshold2 / log2l (x);
    check (x, y);
    y = threshold3 / log2l (x);
    check (x, y);
  }
}

// perform N random tests near overflow threshold
static void
check_near_overflow (int N)
{
  long double threshold1 = 16384.0L;
  long double threshold2 = 16383.0L;
  for (int n = 0; n < N / 2; n++)
  {
    long double x = get_random ();
    x = fabsl (x);
    long double y = threshold1 / log2l (x);
    check (x, y);
    y = threshold2 / log2l (x);
    check (x, y);
  }
}

// perform N random tests for x^y near 1-2^-64, 1-2^-65, 1+2^-64, 1+2^-63
static void
check_near_one (int N)
{
  long double threshold1 = 0x1.fffffffffffffffep-1L; // nextbelow(1) = 1-2^-64
  long double threshold2 = 0x1.0000000000000002p+0L; // nextabove(1) = 1+2^-63
  for (int n = 0; n < N / 4; n++)
  {
    long double x = get_random ();
    x = fabsl (x);
    long double y = threshold1 / log2l (x);
    check (x, y);
    // if x^y is near 1-2^-64, then x^(y/2) is near sqrt(1-2^-64) ~ 1 - 2^-65
    check (x, y * 0.5L);
    y = threshold2 / log2l (x);
    check (x, y);
    // if x^y is near 1+2^-63, then x^(y/2) is near sqrt(1+2^-63) ~ 1 + 2^-64
    check (x, y * 0.5L);
  }
}

// check exact or midpoint 3rd powers in one binade
static void
check_exact_or_midpoint_3 (void)
{
  long double x, y = 3.0L;
  for (x = 1664511.0L; x <= 3329021.0L; x++)
    check (x, y); // check will also try (-x,y)
}

// check exact or midpoint values for y integer
static void
check_exact_or_midpoint (void)
{
  long double zmin = 0x1p-16445L;
  long double zmax = 0x1.fffffffffffffffep+16383L;
  // max_pow[n] is the largest x such that x^n fits in 65 bits
  long double max_pow[] = {0, 0, 6074000999, 3329021, 77935, 8191, 1824, 624, 279, 149, 90, 60, 42, 31, 24, 20, 16, 14, 12, 10, 9, 8, 7, 7, 6, 6, 5, 5, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3};
  // max_m[ey] is the largest m such that m^(2^ey) fits in 64 bits
  long double max_m[] = {0x1.fffffffffffffffep+63L, 0x1.fffffffep+31L,
                         0x1.fffep+15L, 0x1.fep+7L, 0x1.ep+3L, 0x1.8p+1L};
  for (int ey = 5; ey >= 0; ey--)
  {
    int dn = (ey == 0) ? 1 : 2; // for ey > 0, we can restrict to odd n
    // we limit to n >= 4, since n=2 and n=3 take a long time
    int d = 1 << ey; // denominator of y
    for (int n = 41; n >= 4; n -= dn)
    {
      long double y = (long double) n / (long double) d;
      long double xmin = powl (zmin, 1.0L / y);
      long double xmax = powl (zmax, 1.0L / y);
      for (long double m = 3.0L; m <= max_pow[n] && m <= max_m[ey]; m += 2.0L)
      {
        // x = m^d*2^e with m odd and e divisible by d
        long double md = powl (m, d);
        long double tmin = xmin / md;
        long double tmax = xmax / md;
        // we want tmin <= 2^e <= tmax
        int emin, emax;
        frexpl (tmin, &emin); // 2^(emin-1) <= tmin < 2^emin
        frexpl (tmax, &emax); // 2^(emax-1) <= tmax < 2^emax
        // we want emin divisible by d
        while (emin % d) emin++;
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
        for (int e = emin; e <= emax; e += d)
        {
          ref_init();
          ref_fesetround(rnd);
          fesetround(rnd1[rnd]);
          long double x = ldexpl (md, e);
          check (x, y);
        }
      }
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
  ref_fesetround(rnd);
  fesetround(rnd1[rnd]);

  printf ("Checking exact/midpoint 3rd powers\n");
  check_exact_or_midpoint_3 ();

  printf ("Checking exact/midpoint values\n");
  check_exact_or_midpoint ();

  printf ("Checking x=2^k\n");
  check_pow2 ();

#define N 1000000UL /* total number of tests */

  printf ("Checking near overflow threshold\n");
  check_near_overflow (N);

  printf ("Checking near underflow threshold\n");
  check_near_underflow (N);

  printf ("Checking near one\n");
  check_near_one (N);

  printf ("Checking random values\n");

	long int total = 0, fails = 0;

  unsigned int seed = getpid ();
  srand (seed);
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for reduction (+: total,fails)
#endif
	for(uint64_t n = 0; n < N; n++) {
		ref_init();
		ref_fesetround(rnd);
		fesetround(rnd1[rnd]);
		long double x = get_random(), y = get_random();
		int j = check(x, y);
                fails += j;
		total++;
	}
	printf("%ld tests, %ld failure(s)\n", total, fails);

  return 0;
}

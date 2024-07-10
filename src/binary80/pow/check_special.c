#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#include <mpfr.h>
#include <omp.h>
#include <unistd.h>

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

static int
check (long double x, long double y)
{
  long double z, t;
  mpfr_t X, Y, Z;
  mpfr_init2 (X, 64);
  mpfr_init2 (Y, 64);
  mpfr_init2 (Z, 64);
  mpfr_set_ld (X, x, MPFR_RNDN);
  mpfr_set_ld (Y, y, MPFR_RNDN);
  z = cr_powl (x, y);
  t = ref_powl (x, y);
	mpfr_clear (X);
  mpfr_clear (Y);
  mpfr_clear (Z);
	if (!is_equal (z, t))
  {
		if(is_equal(z, -1.)) {
			//printf("GIVEUP x=%La,y=%La (ref = %La)\n", x, y, t);
			return 2;
		} else {
			printf("FAIL x=%La,y=%La ref=%La z=%La\n", x,y,t,z);
#ifdef DO_NOT_ABORT
			return 1;
#else
			exit(1);
#endif
		}
  }
	return 0;
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
  printf ("Checking random values\n");

	long int total = 0, fails = 0, giveups = 0;

#define N 1000000UL /* total number of tests */

  unsigned int seed = getpid ();
  srand (seed);
#pragma omp parallel for reduction (+: total,fails,giveups)
	for(uint64_t n = 0; n < N; n++) {
		ref_init();
		ref_fesetround(rnd);
		fesetround(rnd1[rnd]);
		long double x = get_random(), y = get_random();
		int j = check(x, y);
		if(j == 2) {
			giveups++;
		} else if(j == 1) {
			fails++;
		}
		total++;
	}
	printf("%ld tests, %ld fails, %ld giveups\n", total, fails, giveups);

  return 0;
}

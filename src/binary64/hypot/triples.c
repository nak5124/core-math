/* Generate special cases for hypot testing.

Copyright (c) 2022 Paul Zimmermann, Inria.

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
#include <math.h>
#include <fenv.h>
#include <assert.h>
#include <omp.h>
#include <mpfr.h>

double cr_hypot (double, double);
double ref_hypot (double, double);
extern void ref_init (void);
extern int ref_fesetround (int);

extern int rnd1[];
extern int rnd;
extern int verbose;

char RND[] = "NZUD";

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

static void
doit (double x, double y)
{
  double z1, z2;
  z1 = ref_hypot (x, y);
  fesetround (rnd1[rnd]);
  z2 = cr_hypot (x, y);
  if (asuint64 (z1) != asuint64 (z2)) {
    printf("FAIL x=%la y=%la ref=%la z=%la RND%c\n", x, y, z1, z2, RND[rnd]);
    fflush(stdout);
    exit(1);
  }
}

typedef unsigned __int128 u128;

/* check that x = m * 2^k for 2^52 <= m < 2^53,
   that 2^52 <= y < 2^53, and that z is exactly representable on 54 bits */
static int
valid (u128 x, u128 y, u128 z, int k)
{
  u128 m = x >> k;
  if (x != (m << k))
    return 0;
  if (m < (u128) 0x10000000000000ul || (u128) 0x20000000000000ul <= m)
    return 0;
  if (y < (u128) 0x10000000000000ul || (u128) 0x20000000000000ul <= y)
    return 0;
  assert (z > 0);
  int e = __builtin_ctz (z);
  z = z >> e;
  int ret = z < (u128) 0x40000000000000ul;
  return ret;
}

static u128
gcd (u128 a, u128 b)
{
  while (b != 0)
  {
    u128 r = a % b;
    a = b;
    b = r;
  }
  return a;
}

/* generate all inputs x=j*(p^2-q^2), y=j*(2pq) that satisfy
   2^(52+k) <= x < 2^(53+k), 2^52 <= y < 2^53
   Return the number of generated inputs.
   Performs at most max_loop loops on j.
*/
static unsigned long
generate1 (u128 p, u128 q, int k, unsigned long max_loop)
{
  /* ensure p and q are coprime, otherwise we will get duplicates */
  if (gcd (p, q) != 1)
    return 0;
  unsigned long count = 0;
  assert (p > q);
  u128 x = p * p - q * q;
  u128 y = 2 * p * q;
  u128 z = p * p + q * q;
  u128 xmax = 0x1ffffffffffffful << k; /* (2^53-1)*2^k */
  u128 ymax = 0x1ffffffffffffful;      /* 2^53-1 */
  /* since we want xj = m*2^k >= 2^(52+k), thus j >= 2^(52+k)/x */
  u128 jmin = (u128) 1 << (52+k);
  jmin = (jmin + x - 1) / x;
  for (u128 j = jmin; j <= jmin + max_loop; j++)
  {
    u128 xj = j * x;
    u128 yj = j * y;
    u128 zj = j * z;
    if (xj > xmax)
      break;
    if (yj > ymax)
      break;
    if (valid (xj, yj, zj, k))
    {
      doit (xj, yj);
      count ++;
    }
  }
  return count;
}

/* generate all inputs x=j*(2pq), y=j*(p^2-q^2) that satisfy
   2^(52+k) <= x < 2^(53+k), 2^52 <= y < 2^53
   Return the number of generated inputs. */
static unsigned long
generate2 (u128 p, u128 q, int k, unsigned long max_loop)
{
  /* ensure p and q are coprime, otherwise we will get duplicates */
  if (gcd (p, q) != 1)
    return 0;
  unsigned long count = 0;
  assert (p > q);
  u128 x = 2 * p * q;
  u128 y = p * p - q * q;
  u128 z = p * p + q * q;
  u128 xmax = 0x1ffffffffffffful << k; /* (2^53-1)*2^k */
  u128 ymax = 0x1ffffffffffffful;      /* 2^53-1 */
  /* since we want xj = m*2^k >= 2^(52+k), thus j >= 2^(52+k)/x */
  u128 jmin = (u128) 1 << (52+k);
  jmin = (jmin + x - 1) / x;
  for (u128 j = jmin; j <= jmin + max_loop; j++)
  {
    u128 xj = j * x;
    u128 yj = j * y;
    u128 zj = j * z;
    if (xj > xmax)
      break;
    if (yj > ymax)
      break;
    if (valid (xj, yj, zj, k))
    {
      doit (xj, yj);
      count ++;
    }
  }
  return count;
}

/* Since the check_pythagorean_triples() function tests a huge number of
   values, and will not terminate in reasonable time. you can define
   REDUCE to a large value to test fewer values. */
#define REDUCE 0x8000000000000

/* check all Pythagorean triples x^2 + y^2 = z^2,
   with 2^52 <= y < 2^53, 2^(52+k) <= x < 2^(53+k),
   and z of the form m*2^e with m < 2^54 */
static void
check_pythagorean_triples (int k)
{
  unsigned long count1 = 0, count2 = 0;

  if (verbose)
    fprintf (stderr, "# k=%d\n", k);

  /* Type 1: x = p^2-q^2, y = 2pq, z = p^2+q^2 */
  /* since y = 2pq < 2^53 and q < p, this gives q <= 67108863 */
#pragma omp parallel for schedule(static,1)
  for (u128 q = 1; q <= 67108863; q++)
  {
    u128 p;
    for (p = q + 1; 2 * p * q < 0x20000000000000ul; p += 2 * REDUCE)
      count1 += generate1 (p, q, k, 0xfffffffffffffffful / REDUCE);
  }

  if (verbose)
    fprintf (stderr, "# Type 1: %lu\n", count1);

  /* Type 2: x = 2pq, y = p^2-q^2, z = p^2+q^2, with p even */
  /* since y = p^2-q^2 >= 2*p-1 and y < 2^53, this gives p <= 2^52 */
#pragma omp parallel for
  for (u128 p = 2; p <= 0x10000000000000 / REDUCE; p++)
  {
    /* we want y < 2^53, thus p^2-q^2 < 2^53 thus p^2 - 2^53 < q^2 */
    u128 q, qmin = 1;
    if (p * p < 0x20000000000000ul)
      q = 0;
    else
    {
      q = sqrt (p * p - 0x20000000000000ul);
      while (q * q <= p * p - 0x20000000000000ul)
        q ++;
    }
    if (qmin < q)
      qmin = q;
    if ((p + qmin) % 2 == 0)
      qmin ++; /* ensure p and q have different parities */
    for (q = qmin; q < p; q += 2)
      count2 += generate2 (p, q, k, 0xfffffffffffffffful / REDUCE);
  }
  if (verbose) {
    fprintf (stderr, "# Type 2: %lu\n", count2);
    fprintf (stderr, "# Total: %lu\n", count1 + count2);
  }
}

/* check hypot(x,y) near z, with 27-bit x */
static void
check_bound (mpfr_t z)
{
  mpfr_t x, y, t, zz;
  int ret, rety, ny;
  mpfr_exp_t emax = mpfr_get_emax ();
  mpfr_exp_t emin = mpfr_get_emin ();

  mpfr_init2 (x, 27);
  mpfr_init2 (y, 53);
  mpfr_init2 (t, 2 * mpfr_get_prec (z));
  mpfr_init2 (zz, 2 * mpfr_get_prec (z));
  mpfr_set_emax (mpfr_get_emax_max ());
  mpfr_set_emin (mpfr_get_emin_min ());
  ret = mpfr_sqr (zz, z, MPFR_RNDN);
  assert (ret == 0);
  mpfr_sqrt_ui (t, 2, MPFR_RNDU);
  mpfr_div (x, z, t, MPFR_RNDZ);
  while (1)
  {
    mpfr_hypot (t, x, x, MPFR_RNDZ);
    if (mpfr_cmp (t, z) >= 0)
      break;
    mpfr_nextabove (x);
  }
  /* xmin = ceil(z/sqrt(2)) */
  while (1)
  {
    if (mpfr_cmp (x, z) > 0)
      break;
    mpfr_set_emax (mpfr_get_emax_max ());
    mpfr_set_emin (mpfr_get_emin_min ());
    ret = mpfr_sqr (t, x, MPFR_RNDN);
    assert (ret == 0);
    ret = mpfr_sub (t, zz, t, MPFR_RNDN);
    assert (ret == 0);
    /* we want y^2 near t */
    rety = mpfr_sqrt (y, t, MPFR_RNDZ);
    ny = 2;
    if (rety == 0 || !mpfr_zero_p (y))
    {
      mpfr_nextbelow (y);
      ny ++;
    }
    while (ny-- > 0)
    {
      mpfr_set_emax (emax);
      mpfr_set_emin (emin);
      doit (mpfr_get_d (x, MPFR_RNDN), mpfr_get_d (y, MPFR_RNDN));
      mpfr_nextabove (y);
    }
    mpfr_nextabove (x);
  }
  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (t);
  mpfr_clear (zz);
  /* restore initial values */
  mpfr_set_emax (emax);
  mpfr_set_emin (emin);
}

static void
check_overflow_and_underflow (void)
{
  mpfr_t z;
  mpfr_exp_t emin = mpfr_get_emin ();
  mpfr_set_emin (-1075);
  mpfr_init2 (z, 54);
  if (verbose)
    printf ("Checking near overflow boundary for directed rounding...\n");
  mpfr_set_d (z, 0x1.fffffffffffffp+1023, MPFR_RNDN); /* DBL_MAX */
  check_bound (z);
  if (verbose)
    printf ("Checking near overflow boundary for rounding to nearest...\n");
  mpfr_nextabove (z); /* middle between DBL_MAX and 2^1024 */
  check_bound (z);
  if (verbose)
    printf ("Checking near underflow boundary for directed rounding...\n");
  mpfr_set_ui_2exp (z, 1, -1074, MPFR_RNDN); /* DBL_MIN */
  check_bound (z);
  if (verbose)
    printf ("Checking near underflow boundary for rounding to nearest...\n");
  mpfr_set_ui_2exp (z, 1, -1075, MPFR_RNDN); /* DBL_MIN/2 */
  check_bound (z);
  mpfr_clear (z);
  mpfr_set_emin (emin);
}

void
doloop (int k0, int k1)
{
  ref_init ();
  ref_fesetround (rnd);
  check_overflow_and_underflow ();
  if (verbose)
    printf ("Checking Pythagorean triples...\n");
  for (int k = k0; k <= k1; k++)
    check_pythagorean_triples (k);
}

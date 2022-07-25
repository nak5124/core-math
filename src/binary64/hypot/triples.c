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

double cr_hypot (double, double);
double ref_hypot (double, double);
int ref_fesetround (int);
void ref_init (void);

extern int rnd1[];
extern int rnd;
extern int verbose;

static void
doit (double x, double y)
{
  double z1, z2;
  z1 = ref_hypot (x, y);
  fesetround (rnd1[rnd]);
  z2 = cr_hypot (x, y);
  if (z1 != z2) {
    printf("FAIL x=%la y=%la ref=%la z=%la\n", x, y, z1, z2);
    fflush(stdout);
    exit(1);
  }
}

typedef unsigned __int128 u128;

/* check that x = m * 2*k for 2^52 <= m < 2^53,
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
   Return the number of generated inputs. */
static unsigned long
generate1 (u128 p, u128 q, int k)
{
  /* ensure p and q are coprime, otherwise we will get duplicates */
  if (gcd (p, q) != 1)
    return 0;
  unsigned long count = 0;
  if (p <= q) printf ("p=%lu q=%lu\n", (uint64_t) p, (uint64_t) q);
  assert (p > q);
  u128 x = p * p - q * q;
  u128 y = 2 * p * q;
  u128 z = p * p + q * q;
  u128 xmax = 0x1ffffffffffffful << k; /* (2^53-1)*2^k */
  u128 ymax = 0x1ffffffffffffful;      /* 2^53-1 */
  for (u128 j = 1; ; j++)
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
generate2 (u128 p, u128 q, int k)
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
  for (u128 j = 1; ; j++)
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
#pragma omp parallel for
  for (u128 q = 1; q <= 67108863; q++)
  {
    u128 p;
    for (p = q + 1; 2 * p * q < 0x20000000000000ul; p += 2)
      count1 += generate1 (p, q, k);
  }

  if (verbose)
    fprintf (stderr, "# Type 1: %lu\n", count1);

  /* Type 2: x = 2pq, y = p^2-q^2, z = p^2+q^2, with p even */
  /* since y = p^2-q^2 >= 2*p-1 and y < 2^53, this gives p <= 2^52 */
#pragma omp parallel for
  for (u128 p = 2; p <= 0x10000000000000; p++)
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
      count2 += generate2 (p, q, k);
  }
  if (verbose) {
    fprintf (stderr, "# Type 2: %lu\n", count2);
    fprintf (stderr, "# Total: %lu\n", count1 + count2);
  }
}

void
doloop (int k0, int k1)
{
  ref_init ();
  ref_fesetround (rnd);
  for (int k = k0; k <= k1; k++)
    check_pythagorean_triples (k);
}

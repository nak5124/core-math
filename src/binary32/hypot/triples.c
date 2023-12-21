/* Generate special cases for hypotf testing.

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

float cr_hypotf (float, float);
float ref_hypot (float, float);
int ref_fesetround (int);
void ref_init (void);

extern int rnd1[];
extern int rnd;
extern int verbose;

static inline uint32_t
asuint (float f)
{
  union
  {
    float f;
    uint32_t i;
  } u = {f};
  return u.i;
}

static void
check_aux (float x, float y)
{
  float z1, z2;
  ref_init();
  ref_fesetround(rnd);
  z1 = ref_hypot(x, y);
  fesetround(rnd1[rnd]);
  z2 = cr_hypotf(x, y);
  if (asuint (z1) != asuint (z2)) {
    printf("FAIL x=%a y=%a ref=%a z=%a\n", x, y, z1, z2);
    fflush(stdout);
    exit(1);
  }
}

void
check (float x, float y)
{
  check_aux (x, y);
  check_aux (x, -y);
  check_aux (-x, y);
  check_aux (-x, -y);
  check_aux (y, x);
  check_aux (y, -x);
  check_aux (-y, x);
  check_aux (-y, -x);
}

/* check that x = m * 2*k for 2^23 <= m < 2^24,
   that 2^23 <= y < 2^24, and that z is exactly representable on 25 bits */
static int
valid (uint64_t x, uint64_t y, uint64_t z, int k)
{
  uint64_t m = x >> k;
  if (x != (m << k))
    return 0;
  if (m < 0x800000 || 0x1000000 <= m)
    return 0;
  if (y < 0x800000 || 0x1000000 <= y)
    return 0;
  assert (z > 0);
  int e = __builtin_ctz (z);
  z = z >> e;
  int ret = z < 0x2000000;
  return ret;
}

uint64_t
gcd (uint64_t a, uint64_t b)
{
  while (b != 0)
  {
    uint64_t r = a % b;
    a = b;
    b = r;
  }
  return a;
}

/* generate all inputs x=j*(p^2-q^2), y=j*(2pq) that satisfy
   2^(23+k) <= x < 2^(24+k), 2^23 <= y < 2^24
   Return the number of generated inputs. */
static unsigned long
generate1 (uint64_t p, uint64_t q, int k)
{
  /* ensure p and q are coprime, otherwise we will get duplicates */
  if (gcd (p, q) != 1)
    return 0;
  unsigned long count = 0;
  assert (p > q);
  uint64_t x = p * p - q * q;
  uint64_t y = 2 * p * q;
  uint64_t z = p * p + q * q;
  uint64_t xmax = 0xfffffful << k;
  uint64_t ymax = 0xfffffful;
  for (uint64_t j = 1; ; j++)
  {
    uint64_t xj = j * x;
    uint64_t yj = j * y;
    uint64_t zj = j * z;
    if (xj > xmax)
      break;
    if (yj > ymax)
      break;
    if (valid (xj, yj, zj, k))
    {
      check (xj, yj);
      count ++;
    }
  }
  return count;
}

/* generate all inputs x=j*(2pq), y=j*(p^2-q^2) that satisfy
   2^(23+k) <= x < 2^(24+k), 2^23 <= y < 2^24
   Return the number of generated inputs. */
static unsigned long
generate2 (uint64_t p, uint64_t q, int k)
{
  /* ensure p and q are coprime, otherwise we will get duplicates */
  if (gcd (p, q) != 1)
    return 0;
  unsigned long count = 0;
  assert (p > q);
  uint64_t x = 2 * p * q;
  uint64_t y = p * p - q * q;
  uint64_t z = p * p + q * q;
  uint64_t xmax = 0xfffffful << k;
  uint64_t ymax = 0xfffffful;
  for (uint64_t j = 1; ; j++)
  {
    uint64_t xj = j * x;
    uint64_t yj = j * y;
    uint64_t zj = j * z;
    if (xj > xmax)
      break;
    if (yj > ymax)
      break;
    if (valid (xj, yj, zj, k))
    {
      check (xj, yj);
      count ++;
    }
  }
  return count;
}

/* check all Pythagorean triples x^2 + y^2 = z^2,
   with 2^23 <= y < 2^24, 2^(23+k) <= x < 2^(24+k),
   and z of the form m*2^e with m < 2^25 */
static void
check_pythagorean_triples (int k)
{
  uint64_t p, q;
  unsigned long count1 = 0, count2 = 0;

  if (verbose)
    fprintf (stderr, "# k=%d\n", k);

  /* Type 1: x = p^2-q^2, y = 2pq, z = p^2+q^2 */
  /* since y = 2pq < 2^24 and q < p, this gives q <= 2895 */
#pragma omp parallel for
  for (q = 1; q <= 2895; q++)
    for (p = q + 1; 2 * p * q < 0x1000000ul; p+=2)
      count1 += generate1 (p, q, k);

  if (verbose)
    fprintf (stderr, "# Type 1: %lu\n", count1);

  /* Type 2: x = 2pq, y = p^2-q^2, z = p^2+q^2, with p even */
  /* since y = p^2-q^2 >= 2*p-1 and y < 2^24, this gives p <= 2^23 */
#pragma omp parallel for
  for (p = 2; p <= 0x800000; p++)
  {
    /* we want y < 2^24, thus p^2-q^2 < 2^24 thus p^2 - 2^24 < q^2 */
    uint64_t qmin = 1;
    if (p * p < 0x1000000ul)
      q = 0;
    else
    {
      q = sqrt (p * p - 0x1000000ul);
      while (q * q <= p * p - 0x1000000ul)
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
#pragma omp parallel for
  for (int k = k0; k <= k1; k++)
    check_pythagorean_triples (k);
}

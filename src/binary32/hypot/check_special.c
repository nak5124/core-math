/* Generate special cases for hypotf testing.

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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#include <math.h>
#include <omp.h>
#include <assert.h>

void doloop (int, int);
void check (float, float);
uint64_t gcd (uint64_t, uint64_t);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

/* Check all Pythagorean triples z^2 = x^2 + y^2 with z in the subnormal
   range. We necessarily have x = r^2 - s^2, y = 2*r*s, z = r^2 + s^2
   with gcd(r,s) = 1 and one of r, s even
   (see https://oeis.org/wiki/Pythagorean_triples).
*/
static void
check_triples_subnormal (void)
{
  uint64_t r, s, x, y, z;
  /* the smallest denormal is 2^-149, the smallest normal is 2^-126,
     thus x, y, z are of the form k*2^-149 with k < 2^23. */

  // type I: r is odd
  for (r = 1; r <= 2896; r += 2)
    for (s = 2; s < r; s += 2)
    {
      if (gcd (r, s) == 1)
      {
        x = r * r - s * s;
        y = 2 * r * s;
        z = r * r + s * s;
        if (z > 0x7fffff)
          break;
        // now (x,y,z) is a primitive Pythagorean triple
        for (int n = 1; ; n++)
        {
          uint64_t nn = n * n;
          uint64_t xx = x * nn, yy = y * nn, zz = z * nn;
          if (zz > 0x7fffff)
            break;
          check (ldexpf (xx, -149), ldexpf (yy, -149));
        }
      }
    }

  // type II: r is even
  for (r = 2; r <= 2896; r += 2)
    for (s = 1; s < r; s += 2)
    {
      if (gcd (r, s) == 1)
      {
        x = r * r - s * s;
        y = 2 * r * s;
        z = r * r + s * s;
        if (z > 0x7fffff)
          break;
        // now (x,y,z) is a primitive Pythagorean triple
        for (int n = 1; ; n++)
        {
          uint64_t nn = n * n;
          uint64_t xx = x * nn, yy = y * nn, zz = z * nn;
          if (zz > 0x7fffff)
            break;
          check (ldexpf (xx, -149), ldexpf (yy, -149));
        }
      }
    }
}

typedef struct { uint32_t u, v; } int_pair;
typedef struct { int_pair *l; uint32_t P; int size, alloc; } List_t[1];

static void
list_init (List_t L)
{
  L->l = NULL;
  L->size = L->alloc = 0;
}

static void
list_realloc (List_t L, int size)
{
  L->l = realloc (L->l, size * sizeof (int_pair));
  L->alloc = size;
}

static void
list_add (List_t L, uint32_t u, uint32_t v)
{
  if (L->size == L->alloc)
    list_realloc (L, 2 * L->alloc + 1);
  L->l[L->size].u = u;
  L->l[L->size].v = v;
  L->size++;
}

static void
list_clear (List_t L)
{
  free (L->l);
}

// return 0 <= c < p*q such that c = a mod p and c = b mod q, gcd(p,q) = 1
static uint32_t
crt (uint32_t a, uint32_t b, uint32_t p, uint32_t q)
{
  for (uint32_t c = a; c < p*q; c += p)
    if ((c % q) == b)
      return c;
  assert (0);
}

static void
list_lift (List_t L, uint32_t p)
{
  List_t newL;
  list_init (newL);
  char *squares;
  squares = malloc(p * sizeof (char));
  for (uint32_t u = 0; u < p; u++)
    squares[u] = 0;
  for (uint32_t u = 0; u < p; u++)
  {
    uint32_t t = (u * u) % p;
    squares[t] = 1;
  }
  for (uint32_t u = 0; u < p; u++)
    for (uint32_t v = 0; v < p; v++)
    {
      uint32_t t = (u * u + v * v + p - 1) % p;
      if (squares[t])
      {
        for (int i = 0; i < L->size; i++)
        {
          uint32_t newu = crt (L->l[i].u, u, L->P, p);
          uint32_t newv = crt (L->l[i].v, v, L->P, p);
          // by symmetry we can restrict to u <= v
          if (newu <= newv)
            list_add (newL, newu, newv);
        }
      }
    }
  L->P *= p;
  free (squares);
  free (L->l);
  L->l = newL->l;
  L->size = newL->size;
  L->alloc = newL->alloc;
}

/* Check pairs (x,y) in subnormal range such that x = u*2^-149, y = v*2^-149,
   with u^2 + v^2 = w^2 + 1, u <= v. We force 2 <= u to avoid the trivial
   solutions u=1, v=w. See https://oeis.org/A050796. */
static void
check_triples_subnormal_above (void)
{
  List_t L;
  list_init (L);
  L->P = 1;
  list_add (L, 0, 0);
  /* L is a list of admissible residues (u mod P, v mod P) such that
     u^2 + v^2 - 1 might be a square mod P. */
  list_lift (L, 3);
  list_lift (L, 5);
  list_lift (L, 7);
  list_lift (L, 11);
  list_lift (L, 13);
  // L has 2054282 elements
#pragma omp parallel for schedule(static,1)
  for (int i = 0; i < L->size; i++)
  {
    uint64_t u = L->l[i].u, v = L->l[i].v, p = L->P;
    if (u <= v) // we can skip solutions with v < u by symmetry
    {
      // ensure 2 <= u
      if (u < 2)
      {
        u += p;
        if (v < u)
          v += p;
      }
      for (uint64_t uu = u, vv = v; uu < 0x800000; uu += p, vv += p)
        check (ldexpf (uu, -149), ldexpf (vv, -149));
    }
  }
  list_clear (L);
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

  /* we check triples with exponent difference 0 <= k <= 12 */
  printf ("Checking near-exact subnormal values\n");
  check_triples_subnormal_above ();
  printf ("Checking exact subnormal values\n");
  check_triples_subnormal ();
  printf ("Checking Pythagorean triples\n");
  doloop(0, 12);
  return 0;
}

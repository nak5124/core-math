/* Special checks for exp10m1f.

Copyright (c) 2022-2024 Paul Zimmermann, Inria.

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

float cr_exp10m1f (float);
float ref_exp10m1f (float);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

typedef union { uint32_t n; float x; } union_t;

float
asfloat (uint32_t n)
{
  union_t u;
  u.n = n;
  return u.x;
}

uint32_t
asuint (float f)
{
  union_t u = {.x = f};
  return u.n;
}

// When x is a NaN, returns 1 if x is an sNaN and 0 if it is a qNaN
static inline int issignaling(float x) {
  union_t _x = {.x = x};

  return !(_x.n & (1ull << 22));
}

static inline int is_nan(float x) {
  union_t _x = {.x = x};

  return (((_x.n >> 23) & 0xff) == 0xff) && (_x.n << 9) != 0;
}

/* check for signaling NaN input */
static void
check_signaling_nan (void)
{
  float snan = asfloat (0x7f800001);
  float y = cr_exp10m1f (snan);
  // check that the signaling bit disappeared
  if (!is_nan (y))
  {
    fprintf (stderr, "Error, exp10m1f(snan) should be NaN, got %la=%x\n",
             y, asuint (y));
    exit (1);
  }
  if (issignaling (y))
  {
    fprintf (stderr, "Error, exp10m1f(snan) should be qnan, got snan=%x\n",
             asuint (y));
    exit (1);
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
  ref_init ();
  ref_fesetround (rnd);

  check_signaling_nan ();

  return 0;
}

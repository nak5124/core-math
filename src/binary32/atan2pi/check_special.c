/* Generate special cases for atan2pif testing.

Copyright (c) 2022 Stéphane Glondu and Paul Zimmermann, Inria.

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
#include <string.h>
#include <fenv.h>
#include <math.h>

float cr_atan2pif (float, float);
float ref_atan2pi (float, float);
int ref_fesetround (int);
void ref_init (void);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

static void
check_random ()
{
  float x, y, z, t;
  ref_init ();
  while (1)
  {
    x = ldexp (drand48 (), -149 + (lrand48 () % 277));
    if (lrand48 () & 1)
      x = -x;
    y = ldexp (drand48 (), -149 + (lrand48 () % 277));
    if (lrand48 () & 1)
      y = -y;
    for (rnd = 0; rnd < 4; rnd++)
    {
      ref_fesetround (rnd);
      fesetround (rnd1[rnd]);
      t = ref_atan2pi (y, x);
      z = cr_atan2pif (y, x);
      if (z != t)
      {
        printf("FAIL x=%a y=%a ref=%a z=%a rnd=%d\n", x, y, t, z, rnd);
        exit (1);
      }
    }
  }
}

int
main ()
{
  /* check random values */
  check_random ();
  return 0;
}

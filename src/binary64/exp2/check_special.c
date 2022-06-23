/* Generate special cases for exp2 testing.

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

double cr_exp2 (double);
double ref_exp2 (double);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

static void
check_subnormal (int64_t n, int e)
{
  double x = ldexp ((double) n, e);
  double y1 = ref_exp2 (x);
  fesetround (rnd1[rnd]);
  double y2 = cr_exp2 (x);
  if (y1 != y2)
  {
    printf ("FAIL x=%la ref=%la z=%la\n", x, y1, y2);
    fflush (stdout);
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

  /* check subnormal results */
  /* x0 is the smallest x such that 2^-1075 <= RN(exp2(x)) */
  double x0 = -1075;
  /* x1 is the smallest x such that 2^-1024 <= RN(exp2(x)) */
  double x1 = -1024;
  /* in the [x0,x1) range, floating-point numbers have an integer part
     of 11 bits, thus we multiply by 2^42 to get integers */
  int64_t n0 = ldexp (x0, 42); /* n0 = -4727899999436800 */
  int64_t n1 = ldexp (x1, 42); /* n1 = -4503599627370496 */
#pragma omp parallel for
  for (int64_t n = n0; n < n1; n++)
    check_subnormal (n, -42);
  /* x2 is the smallest x such that 2^-1022 <= RN(exp2(x)) */
  double x2 = -1022;
  /* in the [x1,x2) range, floating-point numbers have an integer part
     of 10 bits, thus we multiply by 2^43 to get integers */
  n1 = ldexp (x1, 43); /* n1 = -9007199254740992, twice as large as above */
  int64_t n2 = ldexp (x2, 43); /* n2 = -8989607068696576 */
#pragma omp parallel for
  for (int64_t n = n1; n < n2; n++)
    check_subnormal (n, -43);

  return 0;
}

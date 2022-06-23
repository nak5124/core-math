/* Check correctness of univariate binary64 function on worst cases.

Copyright (c) 2022 Stéphane Glondu, Inria.

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

#define _POSIX_C_SOURCE 200809L  /* for getline */

#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fenv.h>

#include "function_under_test.h"

double cr_function_under_test (double);
double ref_function_under_test (double);
int ref_fesetround (int);
void ref_init (void);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;

void
doloop(void)
{
  char *buf = NULL;
  size_t buflength = 0;
  ssize_t n;
  double x, z1, z2;
  int count = 0, failures = 0;

  ref_init();
  ref_fesetround(rnd);

  while ((n = getline(&buf, &buflength, stdin)) >= 0) {
    if (n > 0 && buf[0] == '#') continue;
    if (sscanf(buf, "%la", &x) == 1) {
      z1 = ref_function_under_test(x);
      fesetround(rnd1[rnd]);
      z2 = cr_function_under_test(x);
      if (z1 != z2) {
        printf("FAIL x=%la ref=%la z=%la\n", x, z1, z2);
        fflush(stdout);
#ifdef DO_NOT_ABORT
        failures ++;
#else
        exit(1);
#endif
      }
      count++;
    }
  }
  free(buf);

  printf("%d tests passed, %d failure(s)\n", count, failures);
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
      else
        {
          fprintf (stderr, "Error, unknown option %s\n", argv[1]);
          exit (1);
        }
    }

  doloop();
}

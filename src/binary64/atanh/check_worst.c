/* Check correctness of atanh binary64 function on worst cases.

Copyright (c) 2023 Alexei Sibidanov

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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#include <getopt.h>
#include <errno.h>
#include <math.h>
#include <ctype.h>

#include "function_under_test.h"

typedef union {double f; uint64_t u;} b64u64_u;

double cr_function_under_test (double);
double ref_function_under_test (double);
int ref_fesetround (int);
void ref_init (void);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };
int rnd;
FILE *instream;

int transform(double, double*);
int nextarg(double*);
void test();

long parselong(const char *str){
  char *endptr;
  errno = 0;    /* To distinguish success/failure after call */
  long val = strtol(str, &endptr, 0);
  /* Check for various possible errors. */
  if (errno != 0) {
    perror("strtol");
    exit(EXIT_FAILURE);
  }
  if (endptr == str) {
    fprintf(stderr, "No digits were found\n");
    exit(EXIT_FAILURE);
  }
  return val;
}

int main(int argc, char *argv[]){
  static struct option opts[] = {
    { "rndn",       no_argument, 0, 'n'},
    { "rndz",       no_argument, 0, 'z'},
    { "rndu",       no_argument, 0, 'u'},
    { "rndd",       no_argument, 0, 'd'},
    { "help",       no_argument, 0, 'h'},
    {  "rnd", required_argument, 0, 'r'},
    {"input", required_argument, 0, 'i'},
    { "maxf", required_argument, 0, 'm'},
    {      0,                 0, 0,  0 }
  };
  instream = stdin;
  const char *fname = NULL;
  int maxfailures = 10;
  while (1) {
    int ind = 0, c = getopt_long(argc, argv, "nudzhr:i:x:m:", opts, &ind);
    if (c == -1) break;
    switch (c) {
    case 'n': rnd = 0; break;
    case 'z': rnd = 1; break;
    case 'u': rnd = 2; break;
    case 'd': rnd = 3; break;
    case 'm': maxfailures = parselong(optarg); break;
    case 'r':
      rnd = parselong(optarg);
      if(rnd<0||rnd>3) {
	fprintf(stderr, "Rounding mode %d is outside of the range [0,3].\n",rnd);
	exit(EXIT_FAILURE);
      }
      break;
    case 'i':
      fname = optarg;
      break;
    case 'h': /* print a help message (to be written) */  break;
    }
  }
  
  if(fname){
    instream = fopen(fname, "r");
    if (instream == NULL) {
      fprintf(stderr, "Cannot open file %s for reading.\n",fname);
      exit(EXIT_FAILURE);
    }
  }

  test(maxfailures);

  fclose(instream);

  exit(EXIT_SUCCESS);
}

int is_equal (b64u64_u a, b64u64_u b){
  return a.u == b.u;
}

void test(int maxfailures){
  int count = 0, failures = 0;
  ref_init();
  ref_fesetround(rnd);
  fesetround(rnd1[rnd]);
  double x;
  while (nextarg(&x)) {
    b64u64_u zr, zt;
    zr.f = ref_function_under_test(x);
    zt.f = cr_function_under_test(x);
    if (!is_equal (zr, zt)) {
      if(++failures<maxfailures) printf("FAIL x=%a ref=%a z=%a\n", x, zr.f, zt.f);
    }
    ++count;
  }
  printf("%d test arguments, %d succesfully passed and %d failure(s)\n", count, count-failures, failures);
}

int transform(double x, double *out){
  static int first = 1;
  static double px = __builtin_nan("");
  static long k;
  b64u64_u s = {.f = x};
  if (first || px != x) {
    first = 0;
    px = x;
    k = -1;
  }
  if(++k<2){
    s.u |= k<<63;
    *out = s.f;
    return 1;
  } else {
    return 0;
  }
}

int fillbuf(char **buf, size_t *nbuf){
  // fill the buffer
  ssize_t nget;
  while ((nget = getline(buf, nbuf, instream)) != -1) {
    char *pos = *buf, *ncom = strchr(pos, '#');
    if (ncom==pos) continue;
    if (ncom) memset(ncom, 0, *nbuf - (ncom - pos));
     // check that buffer is not empty
    int nonempty = 0;
    for (size_t i = 0, imax = strlen(pos); i<imax; i++) if (!isspace(pos[i])) { nonempty = 1; break;}
    if (nonempty) break;
  }
  if(nget == -1) {
    free(*buf);
    return 0;
  }
  return 1;
}

int nextarg(double *x){
  static int first = 1;
  static double arg = __builtin_nan("");
  if( !first ){
    if (transform(arg, x)) return 1;
  } else {
    first = 0;
  }
  static char *buf = NULL, *pos;
  static size_t nbuf = 0;
  do {
    if(buf==NULL || !*pos){
      int stat = fillbuf(&buf, &nbuf);
      if(stat==0) return 0;
      pos = buf;
    }
    while(*pos){
      char *pos1;
      errno = 0;
      arg = strtod(pos, &pos1);
      if (errno == ERANGE) {transform(arg, x); return 1;} // just zero or infinity
      if (pos1 == pos) {pos++; continue;}
      pos = pos1;
      transform(arg, x); return 1;
    }
  } while(1);
}

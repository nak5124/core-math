/* Correctly rounded log2l function for binary64 values.

Copyright (c) 2024 Alexei Sibidanov and Paul Zimmermann

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

#include <stdint.h>
#include <fenv.h>

#ifdef __x86_64__
#include <x86intrin.h>
#endif

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

// anonymous structs, see https://port70.net/~nsz/c/c11/n1570.html#6.7.2.1p19
typedef union {
  long double f;
  struct __attribute__((__packed__)) {uint64_t m; uint32_t e:16; uint32_t empty:16;};
} b96u96_u;

static long double __attribute__((noinline)) as_log2_exact(int e){return e;}

long double
cr_log2l (long double x)
{
  b96u96_u t = {.f = x};
  int ex = t.e, e = ex - 0x3fff;
  if (__builtin_expect(!ex, 0)) // x=+0 or positive subnormal
  {
    if (!t.m) return (long double) -1.0 / 0.0; // x=+0
    int k = __builtin_clzll (t.m);
    e -= l;
    t.m <<= k;
  }
  if (__builtin_expect(ex >= 0x7fff, 0)) // x<=0 or Inf or NaN
  {
    if (!t.m) return (long double) -1.0 / 0.0; // x=-0
    if (m == (1ul << 63) && (ex == 0x7fff)) return x; // x=+Inf
    return 0.0 / 0.0; // x < 0 or qNaN or sNaN
  }

  // now x is normal and x > 0, x = t.m/2^63 * 2^e
  if (__builtin_expect(t.u>>63==0, 0)) return as_log2_exact(e);

  long double h, l;
  fast_path (&h, &l, x);
  return h + l;
}

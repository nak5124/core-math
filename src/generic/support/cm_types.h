/* Commonly used data types.

Copyright (c) 2024 The CORE-MATH Authors.

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

#ifndef CM_TYPES_H
#define CM_TYPES_H

#include <stdint.h>

#if (defined(__clang__) && __clang_major__ >= 14) || (defined(__GNUC__) && __GNUC__ >= 14)
typedef unsigned _BitInt(128) u128;
typedef _BitInt(128) i128;
#else
typedef unsigned __int128 u128;
typedef __int128 i128;
#endif

typedef union {
  struct {
    u128 r;
    int64_t _ex;
    uint64_t _sgn;
  };
  struct {
    uint64_t lo;
    uint64_t hi;
    int64_t ex;
    uint64_t sgn;
  };
} dint64_t;

typedef union {
  u128 r;
  struct {
    uint64_t l;
    uint64_t h;
  };
} uint128_t;

typedef union {
  double f;
  uint64_t u;
} f64_u;

typedef union {float f; uint32_t u;} b32u32_u;
typedef union {
  double f;
  struct __attribute__((packed)) {uint64_t m:52;uint32_t e:11;uint32_t s:1;};
  uint64_t u;
} b64u64_u;
typedef uint64_t u64;
typedef int64_t i64;
typedef uint16_t ushort;
typedef union {u128 a; u64 b[2];} u128_u;
// typedef union { double f; uint64_t u; } b64u64;

// only the lower 16 bits of e are used
// 1.0 has encoding m=2^63, e=16383
// -1.0 has encoding m=2^63, e=49151
// 2 has encoding m=2^63, e=16384
// +qnan has encoding m=2^63+2^62, e=32767
// -qnan has encoding m=2^63+2^62, e=65535
// +inf has encoding m=2^63, e=32767
// -inf has encoding m=2^63, e=65535
// +snan has encoding m=2^63+2^62-1, e=32767
// -snan has encoding m=2^63+2^62-1, e=65535
typedef union {long double f; struct {uint64_t m; uint16_t e;};} b80u80_t;

typedef union {
  /* Use a little-endian representation.
     FIXME: adapt for big-endian processors. */
  struct {
    u128 rl;
    u128 rh;
    int64_t _ex;
    uint64_t _sgn;
  };
  struct {
    /* for a non-zero qint, hh has its most significant bit set,
       and the significand is in [1, 2):
       x = (-1)^sgn * m * 2^ex
       with m = hh/2^63 + hl/2^127 + lh/2^191 + ll/2^255
    */
    uint64_t ll; /* lower low part */
    uint64_t lh; /* upper low part */
    uint64_t hl; /* lower high part */
    uint64_t hh; /* upper high part */
    int64_t ex;
    uint64_t sgn;
  };
} qint64_t;

// anonymous structs, see https://port70.net/~nsz/c/c11/n1570.html#6.7.2.1p19
typedef union {
  long double f;
  struct __attribute__((__packed__))
  {uint64_t m; uint32_t e:16; uint32_t empty:16;};
} b96u96_u;

#endif  // CM_TYPES_H

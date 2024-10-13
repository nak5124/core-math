#ifndef CM_INTRIN_COMPAT_H
#define CM_INTRIN_COMPAT_H

#ifdef __x86_64__
#  include <x86intrin.h>
#else
#  if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)
#    include <stdint.h>
#    if defined(_MSC_VER)
#      include <arm64intr.h>
#    endif
#  endif
#endif

/************************ Arm _mm_getcsr compatibility ***********************/

// This code emulates the _mm_getcsr SSE intrinsic by reading the FPCR register.
// fegetexceptflag accesses the FPSR register, which seems to be much slower
// than accessing FPCR, so it should be avoided if possible.
// Adapted from sse2neon: https://github.com/DLTcollab/sse2neon
#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)
typedef struct
{
  uint16_t res0;
  uint8_t  res1  : 6;
  uint8_t  bit22 : 1;
  uint8_t  bit23 : 1;
  uint8_t  bit24 : 1;
  uint8_t  res2  : 7;
  uint32_t res3;
} fpcr_bitfield;

inline static unsigned int _mm_getcsr()
{
  union
  {
    fpcr_bitfield field;
    uint64_t value;
  } r;

#if defined(_MSC_VER) && !defined(__clang__)
  r.value = _ReadStatusReg(ARM64_FPCR);
#else
  __asm__ __volatile__("mrs %0, FPCR" : "=r"(r.value));
#endif
  static const unsigned int lut[2][2] = {{0x0000, 0x2000}, {0x4000, 0x6000}};
  return lut[r.field.bit22][r.field.bit23];
}
#endif  // defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)

#if defined(__x86_64__) || defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)
#define CORE_MATH_MM_GETCSR_SUPPORTED
#endif

#endif   // CM_INTRIN_COMPAT_H

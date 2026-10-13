/* Integer-division and error-bound regression checks for tanq.
   Copyright (c) 2026 Chen-Pang He. MIT license; see tanq.c. */
#include "tanq.c"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>

static void require(int ok, const char *what) {
  if (!ok) { fprintf(stderr, "FAIL: %s\n", what); exit(1); }
}

static void set128(mpz_t z, u128 a) {
  u64 words[2] = {(u64)a, (u64)(a >> 64)};
  mpz_import(z, 2, -1, sizeof(u64), 0, 0, words);
}

static void set384(mpz_t z, u384 a) {
  u64 words[6];
  for (int i = 0; i < 3; i++) {
    words[2*i] = a.l[i]; words[2*i+1] = a.l[i] >> 64;
  }
  mpz_import(z, 6, -1, sizeof(u64), 0, 0, words);
}

static u64 state = 1;
static u64 random64(void) {
  state ^= state >> 12; state ^= state << 25; state ^= state >> 27;
  return state * UINT64_C(2685821657736338717);
}
static u128 random128(void) { return ((u128)random64() << 64) | random64(); }

/* The function corpora use normalized frames, so they do not exercise the
   zero/low-limb leading-zero cases or shifts beyond the 384-bit window. */
static void limb_edges(void) {
  const unsigned shifts[] = {0, 1, 63, 64, 65, 127, 128, 129,
                             255, 256, 383, 384, 385};
  mpz_t x, expected, actual;
  mpz_inits(x, expected, actual, NULL);
  for (int bit = -1; bit < 384; bit++) {
    u384 a = {{0, 0, 0}};
    if (bit >= 0) a.l[bit/128] = (u128)1 << (bit%128);
    set384(x, a);
    int bits = mpz_sgn(x) ? (int)mpz_sizeinbase(x, 2) : 0;
    require(lz_384(a) == 384-bits, "384-bit leading zeros");
    if (bit < 128) require(clz128(a.l[0]) == 128-bits, "128-bit leading zeros");
    for (unsigned j = 0; j < sizeof shifts/sizeof *shifts; j++) {
      unsigned shift = shifts[j];
      mpz_fdiv_q_2exp(expected, x, shift);
      set384(actual, shr_384_sat(a, shift));
      require(mpz_cmp(actual, expected) == 0, "384-bit right shift");
      mpz_mul_2exp(expected, x, shift);
      mpz_fdiv_r_2exp(expected, expected, 384);
      set384(actual, shl_384(a, shift));
      require(mpz_cmp(actual, expected) == 0, "384-bit left shift");
    }
  }
  mpz_clears(x, expected, actual, NULL);
  puts("385 leading-zero frames and 10010 shifts passed");
}

static void division(void) {
  mpz_t a, b, q, actual;
  mpz_inits(a, b, q, actual, NULL);
  long lowest = 0, highest = 0;
  for (int i = 0; i < 1000000; i++) {
    u128 aa = random128() | SIGN_MASK, bb = random128() | SIGN_MASK;
    /* Divisors at and beside both ends of the binade, where the reciprocal
       saturates or its fraction wraps, and with extreme limbs. */
    if (i % 13 == 3) bb = SIGN_MASK + random64() % 4;
    if (i % 13 == 4) bb = ~(u128)0 - random64() % 16;
    if (i % 13 == 5) bb |= (u128)UINT64_MAX << 64;
    if (i % 13 == 6) bb |= UINT64_MAX;
    if (i % 13 == 7) bb &= ~(u128)UINT64_MAX;
    if (i % 13 == 8) bb = SIGN_MASK | (random128() >> (1 + random64() % 127));
    if (i % 13 == 9) bb = ~(random128() >> (1 + random64() % 127));
    if (i % 7 == 0) aa = bb;
    if (i % 7 == 1 && bb != ~(u128)0) aa = bb + 1;
    if (i % 7 == 2 && bb != SIGN_MASK) aa = bb - 1;
    if (i % 7 == 3) aa = SIGN_MASK + random64() % 4;
    if (i % 7 == 4) aa = ~(u128)0 - random64() % 4;
    int up;
    u128 got = divide_fast(aa, bb, &up);
    require(got >= SIGN_MASK, "normalized narrow quotient");
    set128(a, aa); set128(b, bb); set128(actual, got);
    mpz_mul_2exp(a, a, 128-up);
    mpz_fdiv_q(q, a, b);
    mpz_sub(actual, actual, q);
    /* Against the floor: above -8 and below 3, or -5 and 2 after halving. */
    long diff = mpz_get_si(actual);
    require(mpz_fits_slong_p(actual) && diff >= (up ? -4 : -7) && diff <= (up ? 1 : 2),
            "128-bit quotient error bound");
    if (diff < lowest) lowest = diff;
    if (diff > highest) highest = diff;
    if (i % 10) continue;
    {
      /* Exercise long carry/borrow chains as well as random limb patterns. */
      u384 x = {{random128(), random128(), aa}};
      u384 y = {{random128(), random128(), bb}};
      if (i % 11 == 0) y.l[0] = y.l[1] = ~(u128)0;
      if (i % 11 == 1) y.l[0] = y.l[1] = 0;
      if (i % 11 == 2) x = y;
      if (i % 11 == 3) {
        y.l[2] |= 1;
        u384 one = {{1, 0, 0}};
        x = sub_384(y, one);
      }
      u384 z = divide_wide(x, y, &up);
      set384(a, x); set384(b, y); set384(actual, z);
      mpz_mul_2exp(a, a, 384-up);
      mpz_fdiv_q(q, a, b);
      require(mpz_cmp(actual, q) == 0, "384-bit quotient");
    }
  }
  mpz_clears(a, b, q, actual, NULL);
  printf("1000000 bounded narrow (floor error %ld to %ld) and 100000 exact wide "
         "quotients passed\n", lowest, highest);
}

/* Samples validate that the implemented pipelines respect the documented
   bounds. The proof is in tanq_error.md; this check is not that proof. */
static void kernel_bounds(void) {
  mpfr_t x, ref, value, error;
  mpfr_inits2(600, x, ref, value, error, (mpfr_ptr)0);
  mpz_t z;
  mpz_init(z);
  double worst = 0;
  for (int i = 0; i < 10000; i++) {
    int e = i % 2 ? (int)(random64() % 16441)-57 : (int)(random64() % 80)-57;
    u128 m = (random128() & MANTISSA_MASK) | IMPLICIT_BIT;
    set128(z, m); mpfr_set_z(x, z, MPFR_RNDN);
    mpfr_mul_2si(x, x, e-112, MPFR_RNDN);
    mpfr_tan(ref, x, MPFR_RNDN);
    int negative = mpfr_signbit(ref);
    mpfr_abs(ref, ref, MPFR_RNDN);
    u128 f[2], sg[2], flip; int ex[2], up;
    if (fast_pair(m, e, f, ex, &flip)) {
      u128 q = divide_fast(f[0], f[1], &up);
      require(!!flip == negative, "fast quotient sign");
      set128(z, q); mpfr_set_z(value, z, MPFR_RNDN);
      mpfr_mul_2si(error, ref, 128-ex[0]+ex[1]-up, MPFR_RNDN);
      mpfr_sub(error, error, value, MPFR_RNDN); mpfr_abs(error, error, MPFR_RNDN);
      require(mpfr_cmp_ui(error, ZIV_GATE) < 0, "fast quotient bound");
      if (mpfr_get_d(error, MPFR_RNDU) > worst) worst = mpfr_get_d(error, MPFR_RNDU);
    }
    u384 w[2];
    if (wide_pair(m, e, w, ex, sg)) {
      u384 q = divide_wide(w[0], w[1], &up);
      require(!!(sg[0] ^ sg[1]) == negative, "wide quotient sign");
      set384(z, q); mpfr_set_z(value, z, MPFR_RNDN);
      mpfr_mul_2si(error, ref, 384-ex[0]+ex[1]-up, MPFR_RNDN);
      mpfr_sub(error, error, value, MPFR_RNDN); mpfr_abs(error, error, MPFR_RNDN);
      mpfr_div_2ui(error, error, 114, MPFR_RNDN);
      require(mpfr_cmp_ui(error, 1) < 0, "wide quotient bound");
    }
  }
  mpz_clear(z); mpfr_clears(x, ref, value, error, (mpfr_ptr)0);
  printf("10000 paired quotient error-bound samples passed; "
         "largest narrow error %.3f of %u units\n", worst, ZIV_GATE);
}

int main(void) {
  limb_edges(); division(); kernel_bounds();
  return 0;
}

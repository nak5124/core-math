/* Integer-division and interval regression checks for tanq.
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

static void get_fixed(mpfr_t r, const uint32_t *a, int n) {
  mpz_t z;
  mpz_init(z);
  mpz_import(z, n, -1, sizeof(uint32_t), 0, 0, a);
  mpfr_set_z(r, z, MPFR_RNDN);
  mpfr_div_2ui(r, r, 32*(n-2), MPFR_RNDN);
  mpz_clear(z);
}

static void intervals(void) {
  for (int n = 10; n <= 34; n = 2*n-2) {
    uint32_t p[n], t[n], s[n], c[n];
    mpfr_t ref, lo, hi, x, err;
    mpfr_inits2(32*n+128, ref, lo, hi, x, err, (mpfr_ptr)0);
    u64 width = big_pio2(p, n);
    get_fixed(lo, p, n);
    mpfr_set_ui(err, width, MPFR_RNDN);
    mpfr_div_2ui(err, err, 32*(n-2), MPFR_RNDN);
    mpfr_add(hi, lo, err, MPFR_RNDN);
    mpfr_const_pi(ref, MPFR_RNDN); mpfr_div_2ui(ref, ref, 1, MPFR_RNDN);
    require(mpfr_cmp(lo, ref) < 0 && mpfr_cmp(ref, hi) < 0, "Machin interval");
    for (int i = 0; i < 1000; i++) {
      memset(t, 0, sizeof t);
      for (int j = 0; j < n-2; j++) t[j] = random64();
      t[n-3] >>= 1; // [0,1/2); plus explicit near-pi/4 samples below
      if (i % 2) t[n-3] = UINT32_C(0xc90fdaa2);
      u64 radius = adaptive_sincos(s, c, t, n);
      get_fixed(x, t, n);
      mpfr_set_ui(err, radius, MPFR_RNDN);
      mpfr_div_2ui(err, err, 32*(n-2), MPFR_RNDN);
      for (int cosine = 0; cosine < 2; cosine++) {
        if (cosine) mpfr_cos(ref, x, MPFR_RNDN); else mpfr_sin(ref, x, MPFR_RNDN);
        get_fixed(lo, cosine ? c : s, n);
        mpfr_sub(ref, ref, lo, MPFR_RNDN); mpfr_abs(ref, ref, MPFR_RNDN);
        require(mpfr_cmp(ref, err) < 0, "Taylor interval");
      }
    }
    mpfr_clears(ref, lo, hi, x, err, (mpfr_ptr)0);
  }
  puts("Machin and 3000 paired Taylor intervals passed");
}

static void rounders(void) {
  const int exponents[] = {-16500, -16495, -16494, -16383, -16382, -1, 0, 1, 16383, 16384};
  const unsigned modes[] = {RM_NEAREST, (unsigned)-1, RM_UP, RM_DOWN};
  const mpfr_rnd_t rounds[] = {MPFR_RNDN, MPFR_RNDZ, MPFR_RNDU, MPFR_RNDD};
  mpz_t a, b;
  mpz_inits(a, b, NULL);
  mpfr_t x, y, z;
  mpfr_inits2(17000, x, y, (mpfr_ptr)0);
  mpfr_init2(z, 113);
  for (unsigned i = 0; i < sizeof exponents/sizeof *exponents; i++)
    for (int offset = -2; offset <= 2; offset++) {
      int e = exponents[i];
      mpz_set_ui(a, 1); mpz_mul_2exp(a, a, 114);
      if (offset < 0) mpz_sub_ui(a, a, -offset); else mpz_add_ui(a, a, offset);
      mpz_set_ui(b, 1);
      if (e < 114) mpz_mul_2exp(b, b, 114-e); else mpz_mul_2exp(a, a, e-114);
      int n = (mpz_sizeinbase(a, 2) + mpz_sizeinbase(b, 2) + 31)/32 + 2;
      uint32_t aa[n], bb[n];
      memset(aa, 0, sizeof aa); memset(bb, 0, sizeof bb);
      mpz_export(aa, NULL, -1, sizeof(uint32_t), 0, 0, a);
      mpz_export(bb, NULL, -1, sizeof(uint32_t), 0, 0, b);
      mpfr_set_emin(-100000); mpfr_set_emax(100000);
      mpfr_set_z(x, a, MPFR_RNDN); mpfr_set_z(y, b, MPFR_RNDN);
      for (int sign = 0; sign < 2; sign++) {
        if (sign) mpfr_neg(x, x, MPFR_RNDN);
        for (int mode = 0; mode < 4; mode++) {
          mpfr_set_emin(-16493); mpfr_set_emax(16384);
          mpfr_clear_flags();
          int inex = mpfr_div(z, x, y, rounds[mode]);
          mpfr_subnormalize(z, inex, rounds[mode]);
          int expected_flags = (mpfr_underflow_p() ? FE_UNDERFLOW : 0)
                             | (mpfr_overflow_p() ? FE_OVERFLOW : 0);
          u128 expected = asu128(mpfr_get_float128(z, MPFR_RNDN));
          int flags;
          u128 actual = round_ratio(aa, bb, n, sign ? SIGN_MASK : 0, modes[mode], &flags);
          require(actual == expected && flags == expected_flags, "rational rounder/flags");
          mpfr_set_emin(-100000); mpfr_set_emax(100000);
        }
      }
    }
  mpz_clears(a, b, NULL); mpfr_clears(x, y, z, (mpfr_ptr)0);
  puts("400 signed rational rounding and range-boundary cases passed");
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
  division(); intervals(); rounders(); kernel_bounds();
  return 0;
}

// extracts e and m such that x = 2^e*m with m odd.
// Assumes x!=0
static inline
void q_extract(int64_t* e, uint64_t* m, long double x) {
	const b80u80_t cvt_x = {.f = x};
	int shift = __builtin_ctzl(cvt_x.m);
	POWL_DPRINTF("shift=%d\n", shift);
	*e = (cvt_x.e&0x7fff) - 16383 - (63 - shift);
	*m = cvt_x.m >> shift;
}

static inline
void q_extract65(int64_t* e, const qint64_t* z) {
	uint64_t l = (z->hl >> 62) | (z->hh << 2), h = z->hh >> 62;
	if(l & 1) {l+=1; if(l == 0) {h++;}}
	// Round to nearest 65 bit number, (ties to away)

	if(l == 0) { // Bit 1 of h has weight z->ex
		*e = z->ex + __builtin_ctzl(h) - 1;
	} else {
		*e = z->ex + __builtin_ctzl(l) - 65;
	}
}

static
bool check_rb(long double x, long double y, const qint64_t* z) {
	uint64_t m;
	int64_t E;
	q_extract(&E, &m, x);
	POWL_DPRINTF("E = %ld\nm = 0x%016lx\n", E, m);


	uint64_t n;
	int64_t F;
	q_extract(&F, &n, y);
	POWL_DPRINTF("F = %ld\nn = 0x%016lx\n", F, n);

	if (m == 1) { // This is exact iff E*n*2^F is in Z, that is 2^-F | E
		if(F >= 0) return true;
		if(F <= -31) return false;
		return !(E & ((1l << (-F)) - 1)); // Lower |F| bits of E must be 0 if exact.
	}

	/* x is not a power of 2 */	
	if(y < 0) return false;
	if(F > 5 || n > 41 || (n * (1 << F) > 41)) return false; // y>41
	if(F < -5) return false; // The only way this was possible was that m == 1

	if(F < 0) {
		if(E & ((1l << (-F)) - 1)) return false; // If the division is not exact

		int64_t off = n*(E >> (-F)); // Cannot overflow given the ranges.
		int64_t G;
		q_extract65(&G, z); // round65(z) = 2^G*k for some odd k
		if(G != off) return false;
	}

	return true;
}

/* Given a not subnormalized such that
	i)  a is hard to round in the current rounding mode
	ii) a approximates a rounding boundary
modifies a in place so that, when rounded to a binary64, the result is
the correctly rounded rounding boundary.
*/ 
static inline
void exactify(qint64_t* a, unsigned rm) {
	if(rm == FE_TONEAREST) {
		// Since we round to even, the correctly rounded result is this one
		if(a->hh&1) {a->hh++;};
		a->hl = a->lh = a->ll = 0;
	} else {
		// In directed rounding modes this is effectively a round to nearest
		if(a->hl>>63) {
			a->hh++;
			if(__builtin_expect(!a->hh, 0)) {
				a->hh = 1ul << 63;
				a->ex++;
			}
		}
		a->hl = a->lh = a->ll = 0;
	}
}

// extracts e and m such that x = 2^e*m with m odd.
// Assumes x!=0
static inline
void q_extract(int64_t* e, uint64_t* m, long double x) {
	const b80u80_t cvt_x = {.f = x};
	int shift = __builtin_ctzl(cvt_x.m);
	int exponent = (cvt_x.e&0x7fff) - 16383;	
	if(exponent == -16383) {exponent++;}

	POWL_DPRINTF("shift=%d\n", shift);
	*e = exponent - (63 - shift);
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
	POWL_DPRINTF("get_hex(R("SAGE_RE" - 2^%ld*%ld))\n",x,E,m);

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
	if(F > 5 || n > 41 || (F >=0 && (n * (1 << F) > 41))) return false; // y>41
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

/* Given a not subnormalized approximating a rounding boundary
modifies a in place so that it exactly represents the rounding boundary.
*/ 
static inline
void exactify(qint64_t* a) {
	if((a->hl>>62) & 1) {
		uint64_t oldhl = a->hl;
		a->hl += (1ul << 63);
		if(a->hl < oldhl) {
			a->hh++;
		}

		if(__builtin_expect(!a->hh, 0)) {
			a->hh = 1ul << 63;
			a->ex++;
		}
	}
	a->hl &= 1ul<<63; a->lh = a->ll = 0;
}

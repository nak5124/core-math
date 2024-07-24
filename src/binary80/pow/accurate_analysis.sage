# based on ../../binary64/pow/qint.sage
def analyze_exp2poly():
	R64  = RealField(64,  rnd="RNDU")
	R128 = RealField(128, rnd="RNDU")
	R192 = RealField(192, rnd="RNDU")
	R256 = RealField(256, rnd="RNDU")
	y = R64(2^-20)
	Q = [R256('0x1.e4cf5158b8e937b4p-28',16).exact_rational(),
 R256('0x1.b5253d395e817dd25174eab9f7d62b5ep-24',16).exact_rational(),
 R256('0x1.62c0223a5c823fd8ffe606e9a039b95p-20',16).exact_rational(),
 R256('0x1.ffcbfc588b0c686b151cc8730e0dc386p-17',16).exact_rational(),
 R256('0x1.430912f86c7876f4b0a8986b2341f3ecc522db683dd0cedb341e93e453d06234p-13',16).exact_rational(),
 R256('0x1.5d87fe78a673110717f69a514bec330576d399be861406bb90dfdff383209fcp-10',16).exact_rational(),
 R256('0x1.3b2ab6fba4e7729ccbbe0b53eeac5072478ea53e63911d92216bb7b020cccd38p-7',16).exact_rational(),
 R256('0x1.c6b08d704a0bf8b33a762bb32bd2dee9eb88e889b40221c4b453bb01af2bb704p-5',16).exact_rational(),
 R256('0x1.ebfbdff82c58ea86f16b06ec9735fcaa3a2751c30ce69d4c2be93bbb135378e6p-3',16).exact_rational(),
 R256('0x1.62e42fefa39ef35793c7673007e5ed5e81e6864ce5316c5b141a2eb7176074b4p-1',16).exact_rational(),
 R256('0x1p+0', 16).exact_rational()]
	assert len(Q) == 11, "len(Q)==11"
	err = dict()
	err[-1] = 2^-261.066
	# mul_qint_11(r, &Q[0], x)
	r = y*Q[0]
	err[0] = (y.ulp()*Q[0])*y^10 # ignored low part of y
	r = R128(r)
	y = R128(y)
	for i in range(1,4+1):
		# add_qint_22(r, &Q[i], r)
		r = Q[i] + r
		erra = R256(r).ulp()*2 # rounding error in the addition
		errb = r.ulp() # error due to the ignored 128 low bits of r_in
		# mul_qint_22(r, r, x)
		r = r*y
		errc = y*r.ulp()
		errd = r*y.ulp() # ignored low parts of the product
		err[i] = (erra+errb+errc+errd)*y^(10 - i)
	r = R256(r)
	y = R256(y)
	for i in range(5, 9+1):
		# add_qint(r, &Q[i], r)
		r = Q[i] + r
		erra = R256(r).ulp()*2 # rounding error in the addition	
		# mul_qint(r, r, x)
		r = r*y
		errb = r.ulp() * 14 # rounding error of the product
		err[i] = (erra+errb)*y^(10 - i)
	# add_qint(r, &Q[10], r)
	r = Q[10] + r
	err[10] = r.ulp()*2
	for i in range(-1, 10 + 1):
		print("err["+str(i)+"]= ", log(err[i])/log(2.))
	tot = add(err[i] for i in [-1..10])
	print("abs. error = ", log(tot)/log(2.))
	print("rel. error = ", log(tot)/log(2.) - 2^-20)

# based on ../../binary64/qint.h
def analyze_logpoly():
	R64  = RealField(64,  rnd="RNDU")
	R128 = RealField(128, rnd="RNDU")
	R192 = RealField(192, rnd="RNDU")
	R256 = RealField(256, rnd="RNDU")
	y = R64(2^-11.999)
	P = [R256('0x1.3703c68998959a7ep-4',16).exact_rational(),
 R256('0x1.484b195526cc202cp-4',16).exact_rational(),
 R256('0x1.5b9ac9b743e9807259bbc2204717d1bp-4',16).exact_rational(),
 R256('0x1.71547652b824e2c3f74362c2ba0a5abp-4',16).exact_rational(),
 R256('0x1.89f3b1694cffdf6e684a821051448ce8p-4',16).exact_rational(),
 R256('0x1.a61762a7aded93f651d3fb8e3ee1c87cp-4',16).exact_rational(),
 R256('0x1.c68f568d317601ce23c4e9633f9f0594p-4',16).exact_rational(),
 R256('0x1.ec709dc3a03fd749fc155223ec56b276p-4',16).exact_rational(),
 R256('0x1.0c9a84994022d285723a2cd20d41ce5df79051fe4d4330949cb230e9d95d849cp-3',16).exact_rational(),
 R256('0x1.2776c50ef9bfe792ca73314d74fb98cc2d7bfaa0f485160883a21b55515b8838p-3',16).exact_rational(),
 R256('0x1.484b13d7c02a8f86a80e36c7d7506f2c4d0df1122ae9cbea9458aceb239285f2p-3',16).exact_rational(),
 R256('0x1.71547652b82fe1777d0ffda0d23a7d11d6ae9a1b51a668b013f4a3d7449f5fdp-3',16).exact_rational(),
 R256('0x1.a61762a7aded93f645c921dc5df9b38219ec861443086f8e129e6c94930fcddp-3',16).exact_rational(),
 R256('0x1.ec709dc3a03fd749fc15522bc2f8a6c27393f1c24e79b11b2d40dcdda82fdb32p-3',16).exact_rational(),
 R256('0x1.2776c50ef9bfe792ca73314d74fb9741788bf77495755d5a7bdc00752aed9344p-2',16).exact_rational(),
 R256('0x1.71547652b82fe1777d0ffda0d23a7d11d6aef551bad2b4b115f70644372ca7d6p-2',16).exact_rational(),
 R256('0x1.ec709dc3a03fd749fc15522bc2f8a6c27393f1c24e6e4641730d9121f43d2fcap-2',16).exact_rational(),
 R256('0x1.71547652b82fe1777d0ffda0d23a7d11d6aef551bad2b4b1164a2cd9a3f406fp-1',16).exact_rational(),
 R256('0x1.71547652b82fe1777d0ffda0d23a7d11d6aef551bad2b4b1164a2cd9a342649p+0',16).exact_rational(),]
	assert len(P) == 19, "len(P)==19"
	err = dict()
	err[-1] = 2^-250.299
	# mul_qint_11(r, x, &P[0])
	r = y*P[0]
	err[0] = (y.ulp()*P[0])*y^18 # ignored low part of y
	r = R128(r)
	y = R128(y)
	for i in range(1,7+1):
		# add_qint_22(r, &P[i], r)
		r = P[i] + r
		erra = R256(r).ulp()*2 # rounding error in the addition
		errb = r.ulp() # error due to the ignored 128 low bits of r_in
		# mul_qint_22(r, r, x)
		r = r*y
		errc = y*r.ulp() # ignored low limbs of r
		errd = 0 # No ignored low part of y !
		err[i] = (erra+errb+errc+errd)*y^(18 - i)
	r = R256(r)
	y = R256(y)
	for i in range(8, 17+1):
		# add_qint(r, &P[i], r)
		r = P[i] + r
		erra = R256(r).ulp()*2 # rounding error in the addition	
		# mul_qint(r, r, x)
		r = r*y
		errb = r.ulp() * 14 # rounding error of the product
		err[i] = (erra+errb)*y^(18 - i)
	# add_qint(r, &P[i], r) (with i = 18)
	r = P[i] + r
	erra = R256(r).ulp()*2 # rounding error in the addition	
	err[18] = erra

	for i in range(-1, 18 + 1):
		print("err["+str(i)+"]= ", log(err[i])/log(2.))
	tot = add(err[i] for i in [-1..10])
	print("abs. error before final product = ", log(tot)/log(2.))
	rel = tot/(log2(1 + 2^-11.999)/2^-11.999)
	rel += 2^-255 * 14 # relative rounding error of the product
	print("total rel. error = ", R64(log(rel)/log(2.)))
	

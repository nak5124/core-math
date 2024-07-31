# based on ../../binary64/pow/qint.sage
# q is the polynomial string output by "sollya accurate_exp2.sollya"
# analyze_exp2poly(q)
# err[-1]=  -261.066000000000
# err[0]=  -290.069724844166
# err[1]=  -328.756371037638
# err[2]=  -304.859494179307
# err[3]=  -281.671700703331
# err[4]=  -257.902813711393
# err[5]=  -360.995874962087
# err[6]=  -337.996874962087
# err[7]=  -315.997874962087
# err[8]=  -293.998874962087
# err[9]=  -271.999874962087
# err[10]=  -254.000000000000
# abs. error =  -253.896573138428
# rel. error =  -253.896572184093
def analyze_exp2poly(q):
	# since all coefficients of Q are positive, we can get an upper
	# bound of each variable by taking y positive and rounding upwards
	R64  = RealField(64,  rnd="RNDU")
	R128 = RealField(128, rnd="RNDU")
	R192 = RealField(192, rnd="RNDU")
	R256 = RealField(256, rnd="RNDU")
	y = R64(2^-19.999)
	R.<x> = RealField(256)[]
	Q = R(q)
	Q = [Q[10-i] for i in [0..10]]
	assert len(Q) == 11, "len(Q)==11"
	for i in [0..10]:
	    assert Q[i]>=0, "Q[i]>=0"
	err = dict()
	err[-1] = 2^-261.066 # polynomial error given by Sollya
	# mul_qint_11(r, &Q[0], x)
	r = y*Q[0]
	err[0] = y.ulp()*Q[0] # ignored low part of y
	err[0] *= y^9         # r has degree 1 now, and degree 10 at the end
	r = R128(r)
	y = R128(y)
	for i in range(1,4+1):
		# check Q[i] is exact representable on 128 bits
		assert R128(Q[i]) == Q[i], "R128(Q[i]) == Q[i]"
		# add_qint_22(r, &Q[i], r)
		r = Q[i] + r
		# Warning: the error of add_qint_22 is bounded by 2 ulp_128(r)
		# not 2 ulp_256(r)
		erra = R128(r).ulp()*2 # rounding error in the addition
                # r has degree i now
		erra *= y^(10-i)
		# mul_qint_22(r, r, x)
		r_in = r
		r = r*y
		errc = r_in.ulp()*y # ignored low part of the product
		errd = r_in*y.ulp() # ignored low part of the product
                # r has degree i+1 now
		errcd = (errc+errd)*y^(10-(i+1))
		err[i] = erra+errcd
	r = R256(r)
	y = R256(y)
	for i in range(5, 9+1):
		# add_qint(r, &Q[i], r)
		r = Q[i] + r
		erra = R256(r).ulp()*2 # rounding error in the addition	
		erra *= y^(10-i)
		# mul_qint(r, r, x)
		r = r*y
		errb = r.ulp() * 14 # rounding error of the product
		errb *= y^(10-(i+1))
		err[i] = erra+errb
	# add_qint(r, &Q[10], r)
	r = Q[10] + r
	err[10] = r.ulp()*2
	for i in range(-1, 10 + 1):
		print("err["+str(i)+"]= ", log(err[i])/log(2.))
	tot = add(err[i] for i in [-1..10])
	print("abs. error = ", log(tot)/log(2.))
	# since |y| <= 2^-19.999, |2^y| >= 2^(-2^-19.999)
	# thus the relative error is bounded by tot/2^(-2^-19.999)
	rel = tot/2^(-y)
	print("rel. error = ", log(rel)/log(2.))

# based on ../../binary64/qint.h
# sage: analyze_logpoly()
# err[-1]=  -250.299000000000
# err[0]=  -294.701160816241
# err[1]=  -334.982999913949
# err[2]=  -322.983999913949
# err[3]=  -310.984999913949
# err[4]=  -298.985999913949
# err[5]=  -286.986999913949
# err[6]=  -274.987999913949
# err[7]=  -262.988999913949
# err[8]=  -376.987536561107
# err[9]=  -364.988536561107
# err[10]=  -352.989536561107
# err[11]=  -340.990536561107
# err[12]=  -328.991536561107
# err[13]=  -316.992536561107
# err[14]=  -303.993536561107
# err[15]=  -291.994536561108
# err[16]=  -279.995536561107
# err[17]=  -266.996536561108
# err[18]=  -255.000000000000
# abs. error before final product =  -250.298781637947
# total rel. error =  -249.998478825644729
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
	

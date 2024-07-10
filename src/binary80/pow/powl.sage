def get_coarsetbl(m = 9, L = 7):
	Rm = RealField(m)
	R  = RealField(106)
	maxmin = -1
	print("static const _Alignas(16)")
	print("lut_t coarse[" + str(2**L) + "] = {")
	for i in range(2**L):
		# x is in range 1 + i/2**L, 1 + (i+1)/2**L
		# we want to minimise |xr - 1| on that interval
		z = 0
		nh = 1 + (i + 1)/2**L
		nl = 1 + i/2**L
		if(1 + i/2**L >= 2**.5):
			z = 1
			nh /= 2
			nl /= 2
		rl = Rm(RealField(m, rnd='RNDD')(1/nh))
		rh = Rm(RealField(m, rnd='RNDU')(1/nl))
		r = rl.nextbelow()
		minr = 20
		minimizer = 0
		while(r <= rh):
			r = r.nextabove()
			valr = max(
				R(abs(r.exact_rational()*nl - 1)),
				R(abs(r.exact_rational()*nh - 1)))
			if(valr < minr):
				minr = valr
				minimizer = r.exact_rational()
		if(rl <= 1 and 1 <= rh):
			minimizer = 1
			minr = max(R(abs(nl - 1)),R(abs(nh - 1)))
		maxmin = max(maxmin, minr)
		print ("{" + get_hex(Rm(minimizer)) + ", "
				+ print_dd(R(-log2(R(minimizer)))) + "," 
				+ str(z) + "},")
	print("};")
	print(get_hex(maxmin))

def get_finetbl(m = 14, L = 7):
	Rm = RealField(m)
	R  = RealField(106)
	maxmin = -1
	print("static const _Alignas(16)")
	print("lut_t fine[" + str(2**L) + "] = {")
	for i in range(2**L):
		# if i < 2**L-1 then x is in range 1 + i/2**(5+L), 1 + (i+1)/2**(5+L)
		# else x is in range 1 - 1/2**7 + (i-2**(L-1))/2**(6+L),
		# 1 - 1/2**7 + (i + 1 - 2**(L - 1))/2**(6+L)
		if(i < 2**(L-1)):
			nh = 1 + (i + 1)/2**(5+L)
			nl = 1 + i/2**(5+L)
		else:
			nh = 1 - 1/2**7 + (i - 2**(L-1) + 1)/2**(6+L)
			nl = 1 - 1/2**7 + (i - 2**(L-1))/2**(6+L)
		rl = Rm(RealField(m, rnd='RNDD')(1/nh))
		rh = Rm(RealField(m, rnd='RNDU')(1/nl))
		r = rl.nextbelow()
		minr = 20
		minimizer = 0
		while(r <= rh):
			r = r.nextabove()
			valr = max(
				R(abs(r.exact_rational()*nl - 1)),
				R(abs(r.exact_rational()*nh - 1))) 
			if(valr < minr):
				minr = valr
				minimizer = r.exact_rational()
		if(rl <= 1 and 1 <= rh):
			minimizer = 1 
			minr = max(R(abs(nl - 1)),R(abs(nh - 1)))
		minlog = R(-log2(minimizer))
		if(R(-log2(R(minimizer))) >= 2**-7):
			z = 1
			minlog = R(-log2(R(minimizer)) - 2**-7)
		else:
			if(R(-log2(R(minimizer))) <= -2**-7):
				z = -1
				minlog = R(-log2(R(minimizer)) + 2**-7)
			else:
				z = 0
		maxmin = max(maxmin, minr)
		print ("{" + get_hex(Rm(minimizer))
				+ ", " + print_dd(minlog) + ", "
				+ str(z) + "}, //" + get_hex(R(nl)) +","+ get_hex(R(nh)) + "(" + hex(i) + ")")
	print("};")
	print(get_hex(maxmin))

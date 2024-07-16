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
		if(minimizer <= 2**(-.5)):
			z = 1
			nh/=2
			nl/=2
			rl*=2
			rh*=2
			minimizer *= 2
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
		#if(R(-log2(R(minimizer))) >= 2**-7):
		#	z = 1
		#	minlog = R(-log2(R(minimizer)) - 2**-7)
		#else:
		#	if(R(-log2(R(minimizer))) <= -2**-7):
		#		z = -1
		#		minlog = R(-log2(R(minimizer)) + 2**-7)
		#	else:
		z = 0
		maxmin = max(maxmin, minr)
		if i // 2**(L-2) == 1:
			print("{0,0,0,0},")
		else:
			print ("{" + get_hex(Rm(minimizer))
					+ ", " + print_dd(minlog) + ", "
					+ str(z) + "}, //" + get_hex(R(nl)) +","+ get_hex(R(nh)) + "(" + hex(i) + ")")
	print("};")
	print(get_hex(maxmin))

def print_bacsel_command(out,y,e,m,t,t0,t1,d,alpha,nthreads):
   y = R64(y)
   command = "./bacsel -rnd_mode all -prec 256 -n 64 -nn 64 -m "
   command += str(m) + " -t " + str(t) + " -t0 " + str(t0)
   command += " -t1 " + str(t1) + " -y " + get_hex(y) + " -e_in " + str(e)
   command += " -d " + str(d) + " -alpha " + str(alpha)
   command += " -nthreads " + str(nthreads)
   if out==None:
      print (command)
   else:
      out.write(command + "\n")

# generates a set of bacsel commands to check worst-cases in the set S
# defined as follows:
# (a) x,y with y integer in [2,41]
#     (41 is the largest integer n such that 3^n fits in 65 bits)
# (b) x,y with y of the form n/2^F with 3 <= n <= 41, n odd, 1 <= F <= 5
#     (F <= 5 because starting from an odd integer of at most 65 bits,
#      you can take at most 5 square roots while keeping an integer)
# t0=2^63+2^62
# check_S(out="/tmp/in",width=2^30,m=55)
# real 2m36.614s, user 27m47.121s on tartine
# *** x,y=0xc.000000000000001p-4,0x2p+0, distance is: 5.421010862e-20 (-m 64)
# check_S(out="/tmp/in",width=2^31,m=64)
# real 3m31.956s, user 28m51.734s
# check_S(out="/tmp/in",width=2^40,m=64)
# real 36m48.652s, user 764m39.762s
def check_S(m=55,t=20,width=2^30,d=2,alpha=2,nthreads=64,out=None):
   if out != None:
      out = open(out,"w")
   R64 = RealField(64)
   t0 = 2^63 + ZZ.random_element(2^63-width)
   t1 = t0 + width
   for n in range(2,41+1):
      y = R64(n)
      print_bacsel_command(out,y,0,m,t,t0,t1,d,alpha,nthreads)
   for F in range(1,5+1):
      for n in range(3,41+1,2):
         y = R64(n/2^F)
         # we need to check 2^F binades
         for e in range(2^F):
            print_bacsel_command(out,y,e,m,t,t0,t1,d,alpha,nthreads)
   if out != None:
      out.close()

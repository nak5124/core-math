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
   command += " -nthreads " + str(nthreads) + " -v"
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
# check_S(out="/tmp/in",width=2^30,m=55) # with t0=2^63+2^62
# *** x,y=0xc.000000000000001p-4,0x2p+0, distance is: 5.421010862e-20 (-m 64)
# check_S(out="/tmp/in",width=2^41,m=64)
# real 19m8.149s, user 1108m52.587s on tartine (64 cores)
# check_S(out="/tmp/in",width=2^42,m=64)
# real 27m38.904s, user 1617m43.970s on tartine (64 cores)
# check_S(out="/tmp/in",width=2^42,m=165,t=22)
# real 23m9.451s, user 1349m14.649s on tartine (64 cores)
# check_S(out="/tmp/in",width=2^42,m=165,t=22,d=3)
# real 13m20.221s, user 759m31.122s on tartine (64 cores)
# check_S(out="/tmp/in",width=2^44,m=165,t=29,d=4) ***
# real 4m25.815s, user 233m15.848s
# check_S(out="/tmp/in",width=2^45,m=165,t=30,d=4)
# real 7m58.439s, user 445m35.797s
# check_S(out="/tmp/in",width=2^44,m=165,t=28,d=5)
# real 6m15.299s, user 343m11.966s
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

def print_xy(out,x,y):
   if abs(x) >= RR(2^16384) or abs(x) < RR(2^-16445):
      return 0
   s = get_hex(x) + "," + get_hex(y)
   if out==None:
      print (s)
   else:
      out.write(s + "\n")
   return 1

# check exact values for x=2^n and y=m/2^k with m odd, k >= 6
# this gives 988254 (x,y) pairs to check
def check_pow2(out=None):
   if out != None:
      out = open(out,"w")
   # if x = 2^n, then x^y = 2^(n*m/2^k) thus n should be multiple of 2^k
   nsols = 0
   for k in range(6,16446):
      # positive n
      for n in range(2^k,16446,2^k):
         assert n%(2^k)==0, "n%(2^k)==0"
         e = n//(2^k)
         # x^y = 2^(e*m)
         for m in range(1,floor(16445/e)+1,2):
            x = RR(2^n)
            y = RR(m/2^k)
            nsols += print_xy(out,x,y)
            nsols += print_xy(out,1/x,-y)
      # negative n
      for n in range(-2^k,-16446,-2^k):
         assert n%(2^k)==0, "n%(2^k)==0"
         e = n//(2^k)
         # x^y = 2^(e*m)
         for m in range(1,floor(16445/abs(e))+1,2):
            x = RR(2^n)
            y = RR(m/2^k)
            nsols += print_xy(out,x,y)
            nsols += print_xy(out,1/x,-y)
   if out != None:
      out.close()
   return nsols

# split binade [2^(e-1),2^e) into k blocks for y
def doit_bacsel(y,e,k,t0=None,t1=None,neg=false):
   if t0==None:
      if neg:
         t0 = -2^64+1
      else:
         t0 = 2^63
   if t1==None:
      if neg:
         t1 = -2^63+1
      else:
         t1 = 2^64
   h = ceil((t1-t0)/k)
   for i in range(k):
      u0 = t0+h*i
      u1 = min(t0+h*(i+1),t1)
      print ("./doit.sh " + str(u0) + " " + str(u1) + " 64 " + str(e) + " 64 24 64 " + y)

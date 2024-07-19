def print_dd(x):
  y = RR(x)
  z = RR(x - y.exact_rational())
  return get_hex(y) + ", " + get_hex(z)

# compute a double-double approximation h+l of -log2(minimizer)
def get_hl(minimizer):
   prec = 200
   h = RR(n(-log2(minimizer),prec))
   H = h.exact_rational()
   l = RR(n(-log2(minimizer) - H,prec))
   L = l.exact_rational()
   err = n(-log2(minimizer) - H - L,prec)
   if h==0:
      assert err==0, "err=0"
   else:
      err = abs(err/h)
   return h, l, err

# see comments in powl_tables.h
# get_coarsetbl()
# 0x1p-7 -107.225895411898 0.496174262004249
# m is the number of bits of r1[i]
# L is the length of the table
def get_coarsetbl(m = 9, L = 7):
	Rm = RealField(m)
	R  = RealField(106)
	T = dict()
	maxmin = -1
	maxerr = 0 # maximal absolute error
	maxhl = 0  # maximal value of |h+l|
	print("static const _Alignas(16)")
	print("lut_t coarse[" + str(2**L) + "] = {")
	for i in range(2**L):
		# x is in range 1 + i/2**L, 1 + (i+1)/2**L
		# we want to minimise |xr - 1| on that interval
		z = 0
		nh = 1 + (i + 1)/2**L # largest x in interval
		nl = 1 + i/2**L       # smallest x in interval
		rl = Rm(RealField(m, rnd='RNDD')(1/nh)) # lower bound for r
		rh = Rm(RealField(m, rnd='RNDU')(1/nl)) # upper bound for r
		minr = 20
		minimizer = 0
		r = rl
		while(r <= rh):
			valr = max(
				R(abs(r.exact_rational()*nl - 1)),
				R(abs(r.exact_rational()*nh - 1)))
			if(valr < minr):
				minr = valr
				minimizer = r.exact_rational()
			r = r.nextabove()
		if(minimizer <= 2**(-.5)): # r <= 1/sqrt(2) thus x >= sqrt(2)
			z = 1 # add a 1/2 prefactor
    # around 1 and 2, fix r to 1 and .5 respectively
		if(rl <= 1 and 1 <= rh):
			minimizer = 1
			minr = max(R(abs(nl - 1)),R(abs(nh - 1)))
		if(rl <= 2 and 2 <= rh):
			minimizer = 1/2
			minr = max(R(abs(nl - 1)),R(abs(nh - 1)))
		maxmin = max(maxmin, minr)
		h, l, err = get_hl (minimizer*2**z) # we have z+h+l ~= -log2(minimizer)
		maxhl = max(maxhl, abs(h+l))
		maxerr = max(maxerr, err)
		T[i] = (h, l)
		print ("{" + get_hex(Rm(minimizer)) + ", "
				+ get_hex(h) + ", " + get_hex(l) + ","
				+ str(z) + "},//" + get_hex(R(minr)))
	print("};")
	print(get_hex(maxmin), log(maxerr)/log(2.), maxhl)
	return T

# m is the number of bits of r2[i]
# L is the length of the table
# get_finetbl()
# 0x1p-12 -107.270397599854 0.0112272554232541
def get_finetbl(m = 13, L = 7):
	Rm = RealField(m)
	R  = RealField(106)
	T = dict()
	maxmin = -1
	print("static const _Alignas(16)")
	print("lut_t fine[" + str(2**L) + "] = {")
	maxerr = 0
	maxhl = 0
	for i in range(2^L):
		# if i < 2^L-1 then x is in range 1 + i/2^(5+L), 1 + (i+1)/2^(5+L)
		# else x is in range 1 - 1/2^7 + (i-2^(L-1))/2^(6+L),
		# 1 - 1/2^7 + (i+1-2^(L-1))/2^(6+L)
		if(i < 2^(L-1)):
			nh = 1 + (i + 1)/2^(5+L)
			nl = 1 + i/2^(5+L)
		else:
			nh = 1 - 1/2^7 + (i - 2^(L-1) + 1)/2^(6+L)
			nl = 1 - 1/2^7 + (i - 2^(L-1))/2^(6+L)
		rl = Rm(RealField(m, rnd='RNDD')(1/nh))
		rh = Rm(RealField(m, rnd='RNDU')(1/nl))
		minr = 20
		minimizer = 0
		r = rl
		while(r <= rh):
			valr = max(
				R(abs(r.exact_rational()*nl - 1)),
				R(abs(r.exact_rational()*nh - 1))) 
			if(valr < minr):
				minr = valr
				minimizer = r.exact_rational()
			r = r.nextabove()
		if(rl <= 1 and 1 <= rh):
			minimizer = 1 
			minr = max(R(abs(nl - 1)),R(abs(nh - 1)))
		minlog = R(-log2(minimizer))
		z = 0
		maxmin = max(maxmin, minr)
		h, l, err = get_hl (minimizer)
		T[i] = (h,l)
		maxerr = max(maxerr, err)
		if i // 2**(L-2) == 1: # unused entries
			print("{0,0,0,0}, // unused")
		else:
			maxhl = max(maxhl, abs(h+l))
			print ("{" + get_hex(Rm(minimizer))
				+ ", " + get_hex(h) + ", " + get_hex(l) + ", "
					+ str(z) + "}, //" + get_hex(R(nl)) +","+ get_hex(R(nh)) + "(" + hex(i) + ")")
	print("};")
	print(get_hex(maxmin), log(maxerr)/log(2.), maxhl)
	return T

def fast_two_sum(a,b):
   assert a==0 or a.ulp()>=b.ulp(), "a==0 or a.ulp()>=b.ulp()"
   rh = a + b
   e = rh - a
   rl = b - e
   return rh, rl

def high_sum(a,bh,bl):
   rh, e = fast_two_sum (a, bh)
   rl = bl + e
   return rh, rl

# analyse_first_high_sum()
# extra_int= 2048.00000000000 r= RNDU i1= 51 err= -105.003076194750
def analyse_first_high_sum():
   T1 = get_coarsetbl()
   maxerr = 0
   for extra_int in range(-16382,16384):
      for r in 'NZUD':
         R = RealField(53,rnd='RND'+r)
         extra_int = R(extra_int)
         for i1 in range(128):
            mlogr1h = R(T1[i1][0])
            mlogr1l = R(T1[i1][1])
            mlogrh, mlogrl = high_sum (extra_int, mlogr1h, mlogr1l)
            a = mlogrh.exact_rational() + mlogrl.exact_rational()
            b = mlogr1h.exact_rational() + mlogr1l.exact_rational()
            b += extra_int.exact_rational()
            err = abs(a-b)/mlogrh
            if err>maxerr:
               maxerr = err
               print ("extra_int=", extra_int, "r=", 'RND'+r, "i1=", i1, "err=", log(err)/log(2.))

def analyse_second_high_sum(extra_int_min=-16382):
   T1 = get_coarsetbl()
   T2 = get_finetbl()
   maxerr = 0
   maxratio = 0
   for extra_int in range(extra_int_min,16384):
      for r in 'NZUD':
         R = RealField(53,rnd='RND'+r)
         extra_int = R(extra_int)
         for i1 in range(128):
            mlogr1h = R(T1[i1][0])
            mlogr1l = R(T1[i1][1])
            mlogrh, mlogrl = high_sum (extra_int, mlogr1h, mlogr1l)
            for i2 in range(128):
               if 32 <= i2 < 64: # unused values
                  continue
               mlogr2h = R(T2[i2][0])
               mlogr2l = R(T2[i2][1])
               b = mlogr2h.exact_rational() + mlogr1l.exact_rational()
               # print ("extra_int=", extra_int, "r=", 'RND'+r, "i1=", i1, "i2=", i2)
               # print ("mlogrh=", get_hex(mlogrh), "mlogr2h=", get_hex(mlogr2h), "mlogr2l=", get_hex(mlogr2l))
               rh, t = fast_two_sum (mlogrh, mlogr2h)
               assert abs(mlogr2l) <= 2^-53*abs(mlogr2h), "abs(mlogr2l) <= 2^-53*abs(mlogr2h)"
               ratio = abs(mlogr2h/rh)
               if ratio>maxratio:
                  maxratio = ratio
                  print ("extra_int=", extra_int, "r=", 'RND'+r, "i1=", i1, "i2=", i2, "ratio=", ratio)
               mlogr2h, mlogr2l = high_sum (mlogrh, mlogr2h, mlogr2l)
               mlogr2l += mlogrl
               a = mlogr2h.exact_rational() + mlogr2l.exact_rational()
               b += mlogr1h.exact_rational() + mlogr1l.exact_rational()
               b += extra_int.exact_rational()
               err = abs(a-b)/mlogr2h
               if err>maxerr:
                  maxerr = err
                  print ("extra_int=", extra_int, "r=", 'RND'+r, "i1=", i1, "i2=", i2, "err=", log(err)/log(2.))

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
# check_S(out="/tmp/in",width=2^45,m=165,t=29,d=4) ***
# real 7m23.300s, user 411m14.889s
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

def get_S():
   S = []
   for n in range(2,41+1):
      y = R64(n)
      S.append(y)
   for F in range(1,5+1):
      for n in range(3,41+1,2):
         y = R64(n/2^F)
         S.append(y)
   return S

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

def a_mul(a,b):
   rh = a*b
   rl = fma(a,b,-rh)
   return rh, rl

def d_mul(ah,al,bh,bl):
   rh, p = a_mul(ah,bh)
   q = fma(al,bh,p)
   rl = fma(ah,bl,q)
   return rh, rl

# return the 'ulp' of the interval x, i.e., max(ulp(t)) for t in x
# this internal routine is used below
def RIFulp(x):
   return max(x.lower().ulp(),x.upper().ulp())

def analyse_polyeval(rnd='RNDN'):
   R = RealField(53,rnd=rnd)   
   err = dict()
   xh = RIF(-2^-11.999,2^-11.999)
   xl = RIF(-2^-64,2^-64)
   ln2invh = RIF(R("0x1.71547652b82fep+0",16))
   ln2invl = RIF(R("0x1.777d0ffda0d24p-56",16))
   # d_mul(&scaleh, &scalel, ln2invh, ln2invl, xh, xl)
   #    a_mul(scaleh, p, ln2invh, xh)
   #    q = fma(ln2invl, xh, p)
   #    scalel = fma(ln2invh, xl, q)
   scaleh = ln2invh*xh
   u = RIFulp(scaleh)
   p = RIF(-u,u) # |p| < ulp(scaleh)
   # a_mul is exact
   q = ln2invl*xh+p
   e1 = RIFulp(q) # maximal error for the 1st fma
   scalel = ln2invh*xl+q
   e2 = RIFulp(scalel) # maximal error for the 2nd fma
   err[0] = e1+e2
   return err[0]

   

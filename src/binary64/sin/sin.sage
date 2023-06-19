# compute table for 1/(2pi)
def computeT():
   x = n(1/(2*pi),2000)
   m = ZZ(floor(x.exact_rational()*2^1280))
   l = m.digits(2^64)
   print ("static const uint64_t T[" + str(len(l)) + "] = {")
   for i in range(len(l)-1,-1,-1):
      print ("   " + hex(l[i]) + ",")
   print ("};")

# given x a 128-bit float, print it
def out_dint(x):
   s, m, e = x.sign_mantissa_exponent()
   while m!=0 and m.nbits()<128:
      m *= 2
      e -= 1
   h, l = divmod(m,2^64)
   print ("  {.hi = " + hex(h) + ", .lo = " + hex(l) + ", .ex = " + str(e+128) + ", .sgn=" + str((1-s)/2) + "},")

# compute table of sin(2*pi*i/2^11) for 0 <= i < 256
def computeS(out=true):
   R = RealField(128)
   if out:
      print ("static const dint64_t S[256] = {")
   S = []
   for i in range(256):
      s = n(sin(2*pi*i/2^11), 512)
      s = R(s)
      if out:
         out_dint(s)
      S.append(s)
   if out:
      print ("};")
   return S

# compute table of cos(2*pi*i/2^11) for 0 <= i < 256
def computeC(out=true):
   R = RealField(128)
   if out:
      print ("static const dint64_t C[256] = {")
   C = []
   for i in range(256):
      s = n(cos(2*pi*i/2^11), 512)
      s = R(s)
      if out:
         out_dint (s)
      C.append(s)
   if out:
      print ("};")
   return C

# implement Algorithm 2.1 from https://www.vinc17.net/research/papers/these.pdf
# finds the smallest distance between y=a*x-b and integers, 0 <= x < N
# return d the smallest distance, and x the corresponding value
def AlgoLefevre(a,b,N):
   x = frac(a)
   y = 1-frac(a)
   d = frac(b)
   u = v = 1
   j = 0
   # invariant: d = j*frac(a) - frac(b) mod 1
   while true:
      if d<x:
         while x<y:
            if u+v>=N:
               break
            y = y-x
            u = u+v
         if u+v>=N:
            break
         x = x-y
         v = v+u
      else:
         d = d-x
         j += v
         while y<x:
            if u+v>=N:
               break
            x = x-y
            v = v+u
         if u+v>=N:
            break
         y = y-x
         u = u+v
      if u+v>=N:
         break
   return d, j

# return x with smallest value of |x/(2pi) cmod 1| for 2^(e-1) <= x < 2^e
def search(e):
   R = RealField(e+200)
   N = 2^52
   b = R(2^(e-1)/(2*pi))
   b = frac(-b)
   if b < 0:
     b = 1+b
   assert b > 0, "b > 0"
   a = R(2^(e-53)/(2*pi))
   a = frac(a)
   assert a > 0, "a > 0"
   oldx = -1
   d, x = AlgoLefevre(a,b,N)
   d, x = AlgoLefevre(a,b+R(d),N)
   x = RR(2^(e-1)+x*2^(e-53))
   return x

def search_all():
   for e in range(1024,1,-1):
      x = search(e)
      print (get_hex(x) + ' # ' + str(e))

# for 2^(e-1) <= x < 2^e
# sin(x) is monotonous between (k-1/2)*pi and (k+1/2)*pi
def doit_bacsel(e):
   x0 = 2^(e-1)
   k0 = ceil(x0/pi+1/2)
   x1 = 2^e
   k1 = floor(x1/pi+1/2)
   t1 = RR(n((k0-1/2)*pi,200))
   if n(t1.exact_rational(),200) < n((k0-1/2)*pi,200):
      t1 = t1.nextabove()
   t1 = ZZ(t1.exact_rational()*2^(53-e))
   print ("./doit.sh 4503599627370496 " + str(t1) + " 53 " + str(e) + " 64 20")
   for k in range(k0+1,k1+1):
      t0 = t1
      t1 = RR(n((k-1/2)*pi,200))
      if n(t1.exact_rational(),200) < n((k-1/2)*pi,200):
         t1 = t1.nextabove()
      t1 = ZZ(t1.exact_rational()*2^(53-e))
      print ("./doit.sh " + str(t0) + " " + str(t1) + " 53 " + str(e) + " 64 20")
   t0 = t1
   print ("./doit.sh " + str(t0) + " 9007199254740992 53 " + str(e) + " 64 20")
   
# return the 'ulp' of the interval x, i.e., max(ulp(t)) for t in x
# this internal routine is used below
def RIFulp(x):
   return max(x.lower().ulp(),x.upper().ulp())

# error analysis of evalPC()
# evalPC()
# err= -125.999957617982
def evalPC(verbose=false):
   RIF128 = RealIntervalField(128)
   X = RIF128(0,2^-11)
   # mul_dint (X2, X, X)
   X2 = X*X
   errX2 = 6*RIFulp(X2)
   PC = [(0x8000000000000000, 0x0, 1, 0),
         (0x9de9e64df22ef2d2, 0x56e26cd9808c1949, 5, 1),
         (0x81e0f840dad61d9a, 0x9980f00630cb655e, 7, 0),
         (0xaae9e3f1e5ffcfe2, 0xa508509534006249, 7, 1),
         (0xf0fa83448dd1e094, 0xe0603ce7044eeba, 6, 0),
	 (0xd368f6f4207cfe49, 0xec63157807ebffa, 5, 1)]
   PC = [2^x[2]*(x[0]/2^64+x[1]/2^128)*(-1)^x[3] for x in PC]
   PC = [RIF128(x) for x in PC]
   # mul_dint_21 (Y, X2, PC+5)
   Y = X2*PC[5]
   err1 = (2*RIFulp(Y)+errX2*PC[5].abs().upper())*X2.abs().upper()^4
   if verbose:
      print ("err1=", log(err1)/log(2.))
   # add_dint (Y, Y, PC+4)
   Y = Y+PC[4]
   err2 = 2*RIFulp(Y)*X2.abs().upper()^4
   if verbose:
      print ("err2=", log(err2)/log(2.))
   # mul_dint (Y, Y, X2)
   Yin = Y
   Y = Y*X2
   err3 = (6*RIFulp(Y)+Yin.abs().upper()*errX2)*X2.abs().upper()^3
   if verbose:
      print ("err3=", log(err3)/log(2.))
   # add_dint (Y, Y, PC+3)
   Y = Y+PC[3]
   err4 = 2*RIFulp(Y)*X2.abs().upper()^3
   if verbose:
      print ("err4=", log(err4)/log(2.))
   # mul_dint (Y, Y, X2)
   Yin = Y
   Y = Y*X2
   err5 = (6*RIFulp(Y)+Yin.abs().upper()*errX2)*X2.abs().upper()^2
   if verbose:
      print ("err5=", log(err5)/log(2.))
   # add_dint (Y, Y, PC+2)
   Y = Y+PC[2]
   err6 = 2*RIFulp(Y)*X2.abs().upper()^2
   if verbose:
      print ("err6=", log(err6)/log(2.))
   # mul_dint (Y, Y, X2)
   Yin = Y
   Y = Y*X2
   err7 = (6*RIFulp(Y)+Yin.abs().upper()*errX2)*X2.abs().upper()
   if verbose:
      print ("err7=", log(err7)/log(2.))
   # add_dint (Y, Y, PC+1)
   Y = Y+PC[1]
   err8 = 2*RIFulp(Y)*X2.abs().upper()
   if verbose:
      print ("err8=", log(err8)/log(2.))
   # mul_dint (Y, Y, X2)
   Yin = Y
   Y = Y*X2
   err9 = 6*RIFulp(Y)+Yin.abs().upper()*errX2
   if verbose:
      print ("err9=", log(err9)/log(2.))
   # add_dint (Y, Y, PC+0)
   Y = Y+PC[0]
   err10 = 2*RIFulp(Y)
   if verbose:
      print ("err10=", log(err10)/log(2.))
   err = err1+err2+err3+err4+err5+err6+err7+err8+err9+err10
   if verbose:
      print ("abs err=", log(err)/log(2.))
   return Y, err

# error analysis of evalPS()
# evalPS()
# (0.01?, 9.1835841043305996799841202077360899805e-41)
def evalPS(verbose=false,xmin=0,xmax=2^-11,rel=false):
   RIF128 = RealIntervalField(128)
   X = RIF128(xmin,xmax)
   # mul_dint (X2, X, X)
   X2 = X*X
   errX2 = 6*RIFulp(X2)
   PS = [(0xc90fdaa22168c234, 0xc4c6628b80dc1cd1, 3, 0),
         (0xa55de7312df295f5, 0x5dc72f712aa57db4, 6, 1),
	 (0xa335e33bad570e92, 0x3f33be0021aa54d2, 7, 0),
	 (0x9969667315ec2d9d, 0xe59d6ab8509a2025, 7, 1),
	 (0xa83c1a43bf1c6485, 0x7d5f8f76fa7d74ed, 6, 0),
	 (0xf16ab2898eae62f9, 0xa7f0339113b8b3c5, 4, 1)]
   PS = [2^x[2]*(x[0]/2^64+x[1]/2^128)*(-1)^x[3] for x in PS]
   PS = [RIF128(x) for x in PS]
   # mul_dint_21 (Y, X2, PS+5)
   Y = X2*PS[5]
   err1 = (2*RIFulp(Y)+errX2*PS[5].abs().upper())*X.abs().upper()^9
   if verbose:
      print ("err1=", log(err1)/log(2.))
   # add_dint (Y, Y, PS+4)
   Y = Y+PS[4]
   err2 = 2*RIFulp(Y)*X.abs().upper()^9
   if verbose:
      print ("err2=", log(err2)/log(2.))
   # mul_dint (Y, Y, X2)
   Yin = Y
   Y = Y*X2
   err3 = (6*RIFulp(Y)+Yin.abs().upper()*errX2)*X.abs().upper()^7
   if verbose:
      print ("err3=", log(err3)/log(2.))
   # add_dint (Y, Y, PS+3)
   Y = Y+PS[3]
   err4 = 2*RIFulp(Y)*X.abs().upper()^7
   if verbose:
      print ("err4=", log(err4)/log(2.))
   # mul_dint (Y, Y, X2)
   Yin = Y
   Y = Y*X2
   err5 = (6*RIFulp(Y)+Yin.abs().upper()*errX2)*X.abs().upper()^5
   if verbose:
      print ("err5=", log(err5)/log(2.))
   # add_dint (Y, Y, PS+2)
   Y = Y+PS[2]
   err6 = 2*RIFulp(Y)*X.abs().upper()^5
   if verbose:
      print ("err6=", log(err6)/log(2.))
   # mul_dint (Y, Y, X2)
   Yin = Y
   Y = Y*X2
   err7 = (6*RIFulp(Y)+Yin.abs().upper()*errX2)*X.abs().upper()^3
   if verbose:
      print ("err7=", log(err7)/log(2.))
   # add_dint (Y, Y, PS+1)
   Y = Y+PS[1]
   err8 = 2*RIFulp(Y)*X.abs().upper()^3
   if verbose:
      print ("err8=", log(err8)/log(2.))
   # mul_dint (Y, Y, X2)
   Yin = Y
   Y = Y*X2
   err9 = (6*RIFulp(Y)+Yin.abs().upper()*errX2)*X.abs().upper()
   if verbose:
      print ("err9=", log(err9)/log(2.))
   # add_dint (Y, Y, PS+0)
   Y = Y+PS[0]
   err10 = 2*RIFulp(Y)*X.abs().upper()
   if verbose:
      print ("err10=", log(err10)/log(2.))
   # mul_dint (Y, Y, X)
   Yin = Y
   Y = Y*X
   err11 = 6*RIFulp(Y)
   if verbose:
      print ("err11=", log(err11)/log(2.))
   err = err1+err2+err3+err4+err5+err6+err7+err8+err9+err10+err11
   if verbose:
      print ("abs err=", log(err)/log(2.))
      print ("Y=", Y.lower(), Y.upper())
   if rel: # convert to relative error
      err = err/Y.abs().lower()
   if verbose:
      print ("rel err=", log(err)/log(2.))
   return Y, err

# relative error for evalPS
# e= -1073 err= 5.9867435683188998247557269010591102197e-38
def evalPSrel():
   maxerr = 0
   for e in range(-11,-1074,-1):
      _, err = evalPS(xmin=2^(e-1),xmax=2^e,rel=true)
      if err>maxerr:
         print ("e=", e, "err=", err)
   return maxerr

# analyze the rounding error when is_sin=1
# analyze_sin_case1()
# i= 255 err= 3.4878355837391092918018785682632564712e-38
# analyze_sin_case1(rel=true)
# i= 0 err= +infinity U= (0.00000000000000000000000000000000000000, 0.0030679615757712824594361751789838895354)
# analyze_sin_case1(rel=true)
# i= 1 err= 1.0823752525137881722193875252693985550e-37 U= (0.0030679423245659126978045966430539338192, 0.0061359039003258701610125534976731766997)
def analyze_sin_case1(rel=false):
   maxerr = 0
   S = computeS(out=false)
   C = computeC(out=false)
   RIF128 = RealIntervalField(128)
   for i in range(1,256):
      U, errU = evalPC()
      V, errV = evalPS()
      # also consider relative error bound of 2^-123.651 from evalPSrel
      errV = min(errV, 2^-123.651*V.abs().upper())
      # mul_dint (U, S+i, U)
      Si = RIF128(S[i])
      Uin = U
      U = Si*U
      err1 = 6*RIFulp(U)+Si.abs().upper()*errU+RIFulp(Si)*Uin.abs().upper()
      # mul_dint (V, C+i, V)
      Ci = RIF128(C[i])
      Vin = V
      V = Ci*V
      err2 = 6*RIFulp(V)+Ci.abs().upper()*errV+RIFulp(Ci)*Vin.abs().upper()
      # add_dint (U, U, V)
      U = U+V
      err3 = 2*RIFulp(U)
      err = err1+err2+err3
      if rel: # convert to relative error
         err = err/U.abs().lower()
      if err>maxerr:
         maxerr = err
         print ("i=", i, "err=", err, "U=", (U.lower(), U.upper()))
      
def repair():
   f = open("torepair","r")
   while true:
      s = f.readline()
      if s=='':
         break
      s = s.split()
      t0 = ZZ(s[0])
      t1 = ZZ(s[1])
      n = ZZ(s[2])
      assert n == 53
      e = ZZ(s[3])
      nthreads = ZZ(s[4])
      assert nthreads == 64
      t = ZZ(s[5])
      assert t == 20
      # find root of sin(x)
      x0 = t0/2^53*2^e
      x1 = t1/2^53*2^e
      k0 = ceil(x0/pi)
      k1 = floor(x1/pi)
      assert k0==k1, "k0==k1"
      tm = round(k0*pi*2^53/2^e)
      print ("./doit.sh " + str(t0) + " " + str(tm-10^7) + " 53 " + str(e) + " 64 20")
      print ("./doit0.sh " + str(tm-10^7) + " " + str(tm+10^7) + " 53 " + str(e) + " 64 20")
      print ("./doit.sh " + str(tm+10^7) + " " + str(t1) + " 53 " + str(e) + " 64 20")
   f.close()

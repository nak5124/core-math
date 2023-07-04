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

def search_naive(e):
   x0 = RR(2^(e-1))
   x1 = RR(2^e)
   k0 = floor(x0/(2*pi))
   k1 = ceil(x1/(2*pi))
   k = max(1,k0)
   mindiff = 1
   best = None
   while k<=k1:
      x = RR(n(k*2*pi,200))
      X = x.exact_rational()
      diff = abs(n(X/(2*pi)-k,200))
      if diff<mindiff:
         best = x
         mindiff = diff
      k += 1
   return best, mindiff

def search_all():
   for e in range(1024,1,-1):
      x = search(e)
      print (get_hex(x) + ' # ' + str(e))

def search2_naive(e):
   x0 = RR(2^(e-1))
   x1 = RR(2^e)
   k0 = floor(x0/pi)
   k1 = ceil(x1/pi)
   k = max(1,k0)
   mindiff = 1
   best = None
   while k<=k1:
      x = RR(n(k*pi,200))
      X = x.exact_rational()
      diff = abs(n(X/pi-k,200))
      if diff<mindiff:
         best = x
         mindiff = diff
      k += 1
   return best, mindiff

# return x with smallest value of |x/pi cmod 1| for 2^(e-1) <= x < 2^e
def search2(e):
   R = RealField(e+200)
   N = 2^52
   b = R(2^(e-1)/pi)
   b = frac(-b)
   if b < 0:
     b = 1+b
   assert b > 0, "b > 0"
   a = R(2^(e-53)/pi)
   a = frac(a)
   assert a > 0, "a > 0"
   oldx = -1
   d, x = AlgoLefevre(a,b,N)
   d, x = AlgoLefevre(a,b+R(d),N)
   x = RR(2^(e-1)+x*2^(e-53))
   return x

def search2_all():
   for e in range(1024,1,-1):
      x = search2(e)
      print (get_hex(x) + ' # ' + str(e))

def out_str(s,out):
   if out==None:
      print(s)
   else:
      out.write(s + '\n')

# for 2^(e-1) <= x < 2^e
# sin(x) is monotonous between (k-1/2)*pi and (k+1/2)*pi
# also avoid roots at k*pi
def doit_bacsel(e,i0=0,i1=infinity,margin=2*10^7,out=None):
   if out != None:
      out = open(out,"w")
   x0 = 2^(e-1)
   k0 = ceil(x0/(pi/2))
   x1 = 2^e
   k1 = floor(x1/(pi/2))
   t1 = RR(n(k0*pi/2,200))
   t1 = ZZ(t1.exact_rational()*2^(53-e))-margin
   i = 0
   if 2^52<t1:
      if i0 <= i < i1:
         out_str ("./doit.sh 4503599627370496 " + str(t1) + " 53 " + str(e) + " 64 20", out)
      i += 1
   else:
      t1 = 2^52
   for k in range(k0+1,k1+2):
      t0 = t1
      t1 = min(t0 + 2*margin,2^53)
      if i0 <= i < i1:
         out_str ("./doit0.sh " + str(t0) + " " + str(t1) + " 53 " + str(e) + " 64 20", out)
      i += 1
      if k==k1+1:
         break
      t0 = t1
      t1 = RR(n(k*pi/2,200))
      t1 = ZZ(t1.exact_rational()*2^(53-e))-margin
      if t0<t1:
         if i0 <= i < i1:
            out_str ("./doit.sh " + str(t0) + " " + str(t1) + " 53 " + str(e) + " 64 20", out)
         i += 1
      else:
         t1 = t0
   t0 = t1
   if t0<2^53:
      if i0 <= i < i1:
         out_str ("./doit.sh " + str(t0) + " 9007199254740992 53 " + str(e) + " 64 20", out)
      i += 1
   if out != None:
      out.close()
   
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

# analyze the rounding error when is_sin=0
# analyze_sin_case2()
# i= 0 err= 6.4652534624868368208105566808251099011e-38 U= (0.99999529380588479549473818088006979171, 1.0000000000000000000000000000000000000)
# analyze_sin_case2(rel=true)
# i= 0 err= 6.4652838893678300901881719744254043831e-38 U= (0.99999529380588479549473818088006979171, 1.0000000000000000000000000000000000000)
def analyze_sin_case2(rel=false):
   maxerr = 0
   S = computeS(out=false)
   C = computeC(out=false)
   RIF128 = RealIntervalField(128)
   for i in range(0,256):
      U, errU = evalPC()
      V, errV = evalPS()
      # also consider relative error bound of 2^-123.651 from evalPSrel
      errV = min(errV, 2^-123.651*V.abs().upper())
      # mul_dint (U, C+i, U)
      Ci = RIF128(C[i])
      Uin = U
      U = Ci*U
      err1 = 6*RIFulp(U)+Ci.abs().upper()*errU+RIFulp(Ci)*Uin.abs().upper()
      # mul_dint (V, S+i, V)
      Si = RIF128(S[i])
      Vin = V
      V = Si*V
      err2 = 6*RIFulp(V)+Si.abs().upper()*errV+RIFulp(Si)*Vin.abs().upper()
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

def a_mul(a,b):
   hi = a*b
   u = RIFulp(hi)
   lo = RIF(-u,u)
   return hi,lo

def s_mul(a,bh,bl):
   hi, lo = a_mul(a,bh)
   lo = a*bl+lo
   return hi, lo

def d_mul(ah,al,bh,bl):
   hi, s = a_mul(ah,bh)
   t = al*bh+s
   lo = ah*bl+t
   return hi,lo, (RIFulp(t)+RIFulp(lo)+RIF(al)*RIF(bl)).abs().upper()

# analyze the absolute error of evalPSfast
# 6.176205449204?e-24
def evalPSfast(verbose=false,xh=None,xl=None):
   P7 = RR("-0x1.33155a7aff959p6",16)
   P5 = RR("0x1.466bc678d8e3fp6", 16)
   P3 = RR("-0x1.4abbce625be53p5", 16)
   P1h = RR("0x1.921fb54442d18p+2", 16)
   P1l = RR("0x1.1a62645458ee1p-52", 16)
   err0 = 2^-77.307
   if xh==None:
      xh = 2^-11+2^-30
      xh = RIF(-xh,xh)
   if xl==None:
      xl = 2^-52.36
      xl = RIF(-xl,xl)
   # a_mul (&uh, &ul, xh, xh)
   uh = xh*xh
   err_uh = RIFulp(uh)
   ul = RIF(-err_uh,err_uh)
   # ul = __builtin_fma (xh + xh, xl, ul)
   ul = (xh+xh)*xl+ul
   # *h = PSfast[4]
   h = RIF(P7)
   # *h = __builtin_fma (*h, uh, PSfast[3])
   hin = h
   h = hin*uh+RIF(P5)
   err1 = (RIFulp(h)+hin*err_uh)*xh.abs().upper()^5
   if verbose:
      print ("err1=", log(err1)/log(2.))
   # *h = __builtin_fma (*h, uh, PSfast[2])
   hin = h
   h = hin*uh+RIF(P3)
   err2 = (RIFulp(h)+hin*err_uh)*xh.abs().upper()^3
   if verbose:
      print ("err2=", log(err2)/log(2.))
   # s_mul (h, l, *h, uh, ul)
   h, l = s_mul (h, uh, ul)
   err3 = RIFulp(l)*xh.abs().upper()
   if verbose:
      print ("err3=", log(err3)/log(2.))
   # fast_two_sum (h, &t, PSfast[0], *h)
   h = P1h + h
   u = RIFulp(h)
   t = RIF(-u,u)
   err4 = h.abs().upper()*2^-105*xh.abs().upper()
   if verbose:
      print ("err4=", log(err4)/log(2.))
   # *l += PSfast[1] + t
   tmp = P1l + t
   l = l + tmp
   err5 = (RIFulp(tmp)+RIFulp(l))*xh.abs().upper()
   if verbose:
      print ("err5=", log(err5)/log(2.))
   # d_mul (h, l, *h, *l, xh, xl)
   h, l, err6 = d_mul(h,l,xh,xl)
   if verbose:
      print ("err6=", log(err6)/log(2.))
   err = err0+err1+err2+err3+err4+err5+err6
   if verbose:
      print ("err=", log(err)/log(2.))
   return err, h, l

# return relative error bound for evalPSfast
# 0x1.7137449123ef6p-26 < x < +Inf
# 0x1.d619ab67d1407p-29 < x/(2pi)
# e= -28 i= 0 err= -77.0995535331724?
def evalPSfast_all():
   maxerr = 0
   for e in [-28..0]:
      # 2^(e-1) <= h < 2^e
      h = RIF(2^(e-1),2^e)
      imin = floor(h.lower()*2^11)
      imax = floor(h.upper()*2^11)
      for i in [imin..imax]:
         xh = RIF(max(0,h.lower()-i*2^-11),min(2^-11,h.upper()-i*2^-11))
         xl = h.upper()*2^-51.64
         err = evalPSfast(xh=xh,xl=xl)
         # convert to relative error
         err = err/(xh+xl).abs().lower()
         if err>maxerr:
            maxerr = err
            print ("e=", e, "i=", i, "err=", log(err)/log(2.))

# analyze the absolute error of	evalPCfast
# evalPCfast(rel=true)
# 8.70261550847?e-22
def evalPCfast(verbose=false,xh=None,xl=None,rel=false):
   P6 = RR("-0x1.55a5c19e443dcp6",16)
   P4 = RR("0x1.03c1f080ad7f9p6", 16)
   P2 = RR("-0x1.3bd3cc9be45dep4", 16)
   P0h = RR("0x1p+0", 16)
   P0l = RR("-0x1.9249c1ep-77", 16)
   err0 = 2^-75.189
   if xh==None:
      xh = 2^-11+2^-30
      xh = RIF(-xh,xh)
   if xl==None:
      xl = 2^-52.36
      xl = RIF(-xl,xl)
   # a_mul (&uh, &ul, xh, xh)
   uh = xh*xh
   err_uh = RIFulp(uh)
   ul = RIF(-err_uh,err_uh)
   # ul = __builtin_fma (xh + xh, xl, ul)
   ul = (xh+xh)*xl+ul
   # *h = PCfast[4]
   h = RIF(P6)
   # *h = __builtin_fma (*h, uh, PCfast[3])
   hin = h
   h = hin*uh+RIF(P4)
   err1 = (RIFulp(h)+hin*err_uh)*xh.abs().upper()^4
   if verbose:
      print ("err1=", log(err1)/log(2.))
   # *h = __builtin_fma (*h, uh, PCfast[2])
   hin = h
   h = hin*uh+RIF(P2)
   err2 = (RIFulp(h)+hin*err_uh)*xh.abs().upper()^2
   if verbose:
      print ("err2=", log(err2)/log(2.))
   # s_mul (h, l, *h, uh, ul)
   h, l = s_mul (h, uh, ul)
   err3 = RIFulp(l)*xh.abs().upper()
   if verbose:
      print ("err3=", log(err3)/log(2.))
   # fast_two_sum (h, &t, PCfast[0], *h)
   h = P0h + h
   u = RIFulp(h)
   t = RIF(-u,u)
   err4 = h.abs().upper()*2^-105
   if verbose:
      print ("err4=", log(err4)/log(2.))
   # *l += PCfast[1] + t
   tmp = P0l + t
   l = l + tmp
   err5 = (RIFulp(tmp)+RIFulp(l))*xh.abs().upper()
   if verbose:
      print ("err5=", log(err5)/log(2.))
   err = err0+err1+err2+err3+err4+err5
   if verbose:
      print ("err=", log(err)/log(2.))
   if rel:
      # convert to relative error
      err = err/(h+l).abs().lower()
   return err, h, l

SC=[
   ("0x0p+0", "0x0p+0", "0x1p+0"),
   ("-0x1.a834cp-45", "0x1.921f8bec23b6p-9", "0x1.ffff621621d0ap-1"),
   ("-0x1.cef8p-44", "0x1.921f0fe5ba394p-8", "0x1.fffd8858e8ab6p-1"),
   ("-0x1.a8d5p-46", "0x1.2d96b0e4f495dp-7", "0x1.fffa72c978c5bp-1"),
   ("-0x1.e8823p-41", "0x1.921d1fcaed2e2p-7", "0x1.fff62169b9536p-1"),
   ("0x1.e9d64p-42", "0x1.f6a296ad1a28ap-7", "0x1.fff0943c53a57p-1"),
   ("-0x1.db31cp-43", "0x1.2d936bbdd3a6p-6", "0x1.ffe9cb44b527dp-1"),
   ("-0x1.4030ep-42", "0x1.5fd4d21f2d6cbp-6", "0x1.ffe1c6870ccd1p-1"),
   ("0x1.ec708p-43", "0x1.92155f7a97112p-6", "0x1.ffd886084cbddp-1"),
   ("0x1.ccef8p-40", "0x1.c454f4d127745p-6", "0x1.ffce09ce29c7ap-1"),
   ("0x1.31659cp-38", "0x1.f69373249ae67p-6", "0x1.ffc251df1b68ap-1"),
   ("-0x1.0a03cp-38", "0x1.14685db0e8dc1p-5", "0x1.ffb55e42619e1p-1"),
   ("-0x1.6d229p-39", "0x1.2d8657570832dp-5", "0x1.ffa72efff0c79p-1"),
   ("0x1.8efe09p-36", "0x1.46a3971317fd8p-5", "0x1.ff97c4207f829p-1"),
   ("0x1.7d1b7p-40", "0x1.5fc00d2a37dfbp-5", "0x1.ff871dadb7505p-1"),
   ("-0x1.6fb75p-40", "0x1.78dbaa5753e8fp-5", "0x1.ff753bb1b9eadp-1"),
   ("-0x1.860bdp-39", "0x1.91f65f0e798f5p-5", "0x1.ff621e3798b8ep-1"),
   ("-0x1.1f85ap-39", "0x1.ab101bd4352b3p-5", "0x1.ff4dc54b1d65ep-1"),
   ("-0x1.950eep-40", "0x1.c428d12acfd89p-5", "0x1.ff3830f8d68ebp-1"),
   ("0x1.19ac5p-38", "0x1.dd406f9b7c52cp-5", "0x1.ff21614e0fe5fp-1"),
   ("0x1.80b7p-40", "0x1.f656e7a0afa3fp-5", "0x1.ff095658e5f26p-1"),
   ("0x1.bb805p-39", "0x1.07b614e5bea07p-4", "0x1.fef01028234b7p-1"),
   ("-0x1.04b098p-38", "0x1.1440134bd80c4p-4", "0x1.fed58ecb6abp-1"),
   ("0x1.d4324p-40", "0x1.20c9674f8babep-4", "0x1.feb9d25302722p-1"),
   ("-0x1.272734p-36", "0x1.2d520925a8254p-4", "0x1.fe9cdad029914p-1"),
   ("-0x1.cf2eap-40", "0x1.39d9f12ba4ce5p-4", "0x1.fe7ea8548494p-1"),
   ("-0x1.42bfap-39", "0x1.46611791755b4p-4", "0x1.fe5f3af2e61a6p-1"),
   ("0x1.c7aa8p-42", "0x1.52e774a501658p-4", "0x1.fe3e92be9d11fp-1"),
   ("0x1.8f706p-39", "0x1.5f6d00aae2d16p-4", "0x1.fe1cafcbd2534p-1"),
   ("0x1.0f9a6p-40", "0x1.6bf1b3e8054f1p-4", "0x1.fdf9922f72013p-1"),
   ("0x1.6e2e4cp-36", "0x1.787586aec8bd5p-4", "0x1.fdd539ff04d69p-1"),
   ("-0x1.7f29p-42", "0x1.84f8712bed97bp-4", "0x1.fdafa75145ab1p-1"),
   ("-0x1.79566p-40", "0x1.917a6bc207cb5p-4", "0x1.fd88da3d14232p-1"),
   ("-0x1.481a8p-38", "0x1.9dfb6eb049c72p-4", "0x1.fd60d2da7c4ccp-1"),
   ("0x1.c3dp-43", "0x1.aa7b7244abcfp-4", "0x1.fd379142206e6p-1"),
   ("0x1.a74f4p-38", "0x1.b6fa6ec6247f4p-4", "0x1.fd0d158d7d201p-1"),
   ("-0x1.014fp-41", "0x1.c3785c79b9f66p-4", "0x1.fce15fd6db19ep-1"),
   ("-0x1.601ap-42", "0x1.cff533b2e583cp-4", "0x1.fcb4703914b29p-1"),
   ("0x1.1f77bp-38", "0x1.dc70ecbcaa797p-4", "0x1.fc8646cfe4e15p-1"),
   ("-0x1.cbcdp-40", "0x1.e8eb7fdd975e3p-4", "0x1.fc56e3b7dc611p-1"),
   ("-0x1.2dfacp-36", "0x1.f564e5633c113p-4", "0x1.fc26470e37059p-1"),
   ("-0x1.0548ap-39", "0x1.00ee8ad695ba2p-3", "0x1.fbf470f0ac105p-1"),
   ("-0x1.f4fedp-38", "0x1.072a047a21f9dp-3", "0x1.fbc1617e50bc5p-1"),
   ("0x1.75bf8p-40", "0x1.0d64dbcb6f37cp-3", "0x1.fb8d18d66871ap-1"),
   ("-0x1.c8548p-40", "0x1.139f0ced568e7p-3", "0x1.fb5797196077dp-1"),
   ("-0x1.e70ccp-39", "0x1.19d8940b24dccp-3", "0x1.fb20dc6823e96p-1"),
   ("-0x1.67e3bp-38", "0x1.20116d4dafe44p-3", "0x1.fae8e8e476ec3p-1"),
   ("-0x1.b2992p-39", "0x1.264994df2a5b7p-3", "0x1.faafbcb0d5ff4p-1"),
   ("0x1.1bd37p-38", "0x1.2c8106e9c2949p-3", "0x1.fa7557f082241p-1"),
   ("-0x1.fde3bp-37", "0x1.32b7bf9139844p-3", "0x1.fa39bac7bf75fp-1"),
   ("-0x1.08e9ep-37", "0x1.38edbb0b3d941p-3", "0x1.f9fce55aeb11dp-1"),
   ("-0x1.230a2p-39", "0x1.3f22f57d43a42p-3", "0x1.f9bed7cfc2566p-1"),
   ("-0x1.d28a8p-40", "0x1.45576b1239748p-3", "0x1.f97f924c943dp-1"),
   ("-0x1.13ecap-38", "0x1.4b8b17f6c9ce7p-3", "0x1.f93f14f86385cp-1"),
   ("0x1.7f018p-41", "0x1.51bdf859a1755p-3", "0x1.f8fd5ffae290dp-1"),
   ("-0x1.3e262p-38", "0x1.57f00864566acp-3", "0x1.f8ba737cbf352p-1"),
   ("-0x1.61388p-41", "0x1.5e21444891d1bp-3", "0x1.f8764fa71636p-1"),
   ("0x1.1ecd5p-38", "0x1.6451a832b6025p-3", "0x1.f830f4a402949p-1"),
   ("0x1.f5b61p-37", "0x1.6a8130526c4f9p-3", "0x1.f7ea629e40f74p-1"),
   ("-0x1.e8064p-39", "0x1.70afd8cfcfcbfp-3", "0x1.f7a299c1abc2bp-1"),
   ("0x1.ef77cp-40", "0x1.76dd9de56b974p-3", "0x1.f7599a3a0d93dp-1"),
   ("0x1.ada62p-39", "0x1.7d0a7bbdd2789p-3", "0x1.f70f6434b0126p-1"),
   ("-0x1.d726p-43", "0x1.83366e89baf17p-3", "0x1.f6c3f7df5c476p-1"),
   ("0x1.1b401p-37", "0x1.8961727df625p-3", "0x1.f67755686e715p-1"),
   ("-0x1.250cp-40", "0x1.8f8b83c661f19p-3", "0x1.f6297cff78997p-1"),
   ("0x1.d4e02p-37", "0x1.95b49e9e34991p-3", "0x1.f5da6ed4120c2p-1"),
   ("-0x1.0db04p-38", "0x1.9bdcbf2cf4ba2p-3", "0x1.f58a2b17948e7p-1"),
   ("0x1.474318p-34", "0x1.a203e1c13d097p-3", "0x1.f538b1fa20f7p-1"),
   ("-0x1.9634ap-38", "0x1.a82a0259c8277p-3", "0x1.f4e603b0c377ap-1"),
   ("0x1.9ef77p-36", "0x1.ae4f1d643628cp-3", "0x1.f492206b8640dp-1"),
   ("-0x1.ec44p-41", "0x1.b4732ef3a73a6p-3", "0x1.f43d085ffbc0fp-1"),
   ("-0x1.c09888p-36", "0x1.ba963349b5dap-3", "0x1.f3e6bbc207ea3p-1"),
   ("-0x1.2bbe1p-37", "0x1.c0b826a619915p-3", "0x1.f38f3ac66822dp-1"),
   ("0x1.2274p-38", "0x1.c6d90536b5ba6p-3", "0x1.f33685a39e448p-1"),
   ("0x1.3fc51p-36", "0x1.ccf8cb34fdf74p-3", "0x1.f2dc9c9051226p-1"),
   ("-0x1.05094p-39", "0x1.d31774d268105p-3", "0x1.f2817fc466752p-1"),
   ("0x1.64a38p-39", "0x1.d934fe54dc742p-3", "0x1.f2252f775b969p-1"),
   ("0x1.51da3p-37", "0x1.df5163f2148edp-3", "0x1.f1c7abe265636p-1"),
   ("0x1.8aabp-39", "0x1.e56ca1e198337p-3", "0x1.f168f53f68d6cp-1"),
   ("-0x1.2c34ep-38", "0x1.eb86b461f9509p-3", "0x1.f1090bc8a71b9p-1"),
   ("0x1.88c1ep-38", "0x1.f19f97b3412bap-3", "0x1.f0a7efb910508p-1"),
   ("0x1.11cf1p-37", "0x1.f7b7480d7462ep-3", "0x1.f045a14cdcc68p-1"),
   ("-0x1.1ea75p-36", "0x1.fdcdc1aa96acp-3", "0x1.efe220c0f169cp-1"),
   ("-0x1.7bcc98p-36", "0x1.01f180695e83ap-2", "0x1.ef7d6e52155fbp-1"),
   ("-0x1.590bp-41", "0x1.04fb80e36f7a2p-2", "0x1.ef178a3e4964bp-1"),
   ("-0x1.2b2cp-41", "0x1.0804e05ea8318p-2", "0x1.eeb074c50c38fp-1"),
   ("-0x1.fb9f4p-38", "0x1.0b0d9cfcfd474p-2", "0x1.ee482e25c3da2p-1"),
   ("0x1.042ccp-37", "0x1.0e15b4e239b7cp-2", "0x1.eddeb6a05d726p-1"),
   ("0x1.84f19p-37", "0x1.111d262c45d07p-2", "0x1.ed740e765bd98p-1"),
   ("-0x1.22a4ep-38", "0x1.1423eefbfb4f4p-2", "0x1.ed0835e9a8644p-1"),
   ("-0x1.d144p-38", "0x1.172a0d76b54d9p-2", "0x1.ec9b2d3c54ep-1"),
   ("0x1.529b1p-37", "0x1.1a2f7fbf8ec8ep-2", "0x1.ec2cf4b18ac69p-1"),
   ("0x1.7ed7ap-38", "0x1.1d3443f55e187p-2", "0x1.ebbd8c8ddbc78p-1"),
   ("-0x1.30aeb8p-36", "0x1.2038583ba73d1p-2", "0x1.eb4cf515fbdbap-1"),
   ("-0x1.f9486p-38", "0x1.233bbabb7d7c4p-2", "0x1.eadb2e8e96c05p-1"),
   ("0x1.6b72p-39", "0x1.263e699599a5fp-2", "0x1.ea68393e5b3f4p-1"),
   ("-0x1.20e4cp-37", "0x1.294062ec80dp-2", "0x1.e9f4156c83cc5p-1"),
   ("-0x1.4ec28p-37", "0x1.2c41a4e858f51p-2", "0x1.e97ec3603d3efp-1"),
   ("-0x1.64f65p-37", "0x1.2f422dadb4708p-2", "0x1.e90843620902cp-1"),
   ("0x1.5a13ccp-35", "0x1.3241fb6799231p-2", "0x1.e89095ba336c9p-1"),
   ("-0x1.c3792p-38", "0x1.35410c2d6f116p-2", "0x1.e817bab4e7d66p-1"),
   ("-0x1.d2d84p-39", "0x1.383f5e34e41e1p-2", "0x1.e79db29a5f5f6p-1"),
   ("-0x1.38398p-40", "0x1.3b3cefa024218p-2", "0x1.e7227db6ae2c2p-1"),
   ("-0x1.492d1p-36", "0x1.3e39be9500afbp-2", "0x1.e6a61c5625928p-1"),
   ("-0x1.e62bcp-39", "0x1.4135c9411bbbep-2", "0x1.e6288ec49d09fp-1"),
   ("-0x1.b99e9p-37", "0x1.44310dc74a6d9p-2", "0x1.e5a9d5507d64dp-1"),
   ("0x1.5c2b08p-36", "0x1.472b8a5777419p-2", "0x1.e529f046d2a19p-1"),
   ("-0x1.73f78p-39", "0x1.4a253d11730c8p-2", "0x1.e4a8dff828abfp-1"),
   ("-0x1.5020cp-39", "0x1.4d1e2427500ep-2", "0x1.e426a4b2c6d41p-1"),
   ("-0x1.14bcp-41", "0x1.50163dc18a2f9p-2", "0x1.e3a33ec75f23p-1"),
   ("-0x1.12634p-37", "0x1.530d880a28693p-2", "0x1.e31eae87308fbp-1"),
   ("-0x1.1026cp-38", "0x1.5604012ee1bedp-2", "0x1.e298f443a370bp-1"),
   ("-0x1.00ba4p-39", "0x1.58f9a75a8287ap-2", "0x1.e212104f70ecp-1"),
   ("0x1.c6843p-36", "0x1.5bee78bc7ab57p-2", "0x1.e18a02fd4d22cp-1"),
   ("0x1.21f708p-36", "0x1.5ee2737b964fep-2", "0x1.e100cca24a015p-1"),
   ("0x1.b06eap-37", "0x1.61d595c9cad64p-2", "0x1.e0766d924647ap-1"),
   ("-0x1.b333p-40", "0x1.64c7ddd3ca7p-2", "0x1.dfeae622e3542p-1"),
   ("0x1.65b64cp-35", "0x1.67b949cef2662p-2", "0x1.df5e36a8f4f63p-1"),
   ("0x1.b12cp-43", "0x1.6aa9d7dc7cda1p-2", "0x1.ded05f7de38cap-1"),
   ("0x1.1480bp-36", "0x1.6d99863a367f9p-2", "0x1.de4160f68b4f4p-1"),
   ("0x1.93ff9p-37", "0x1.70885310cc63bp-2", "0x1.ddb13b6c930aep-1"),
   ("-0x1.40e94p-38", "0x1.73763c91eb992p-2", "0x1.dd1fef38bff12p-1"),
   ("0x1.c2dbep-37", "0x1.766340f38b25fp-2", "0x1.dc8d7cb3cf6a8p-1"),
   ("0x1.f1de8p-40", "0x1.794f5e616b6b5p-2", "0x1.dbf9e4394e586p-1"),
   ("0x1.0ce78p-40", "0x1.7c3a9311f5519p-2", "0x1.db6526238522fp-1"),
   ("0x1.39238p-37", "0x1.7f24dd3818314p-2", "0x1.dacf42ce3aa8dp-1"),
   ("-0x1.163d6p-38", "0x1.820e3b0485789p-2", "0x1.da383a967d314p-1"),
   ("0x1.0ebcep-37", "0x1.84f6aaaffdb71p-2", "0x1.d9a00dd88b71ep-1"),
   ("-0x1.30f98p-35", "0x1.87de2a6775698p-2", "0x1.d906bcf3e027cp-1"),
   ("0x1.7fcb88p-35", "0x1.8ac4b871b75bp-2", "0x1.d86c484371db1p-1"),
   ("-0x1.3cd02p-37", "0x1.8daa52eba4ff7p-2", "0x1.d7d0b02bbf203p-1"),
   ("-0x1.7b168p-38", "0x1.908ef81e6ebbbp-2", "0x1.d733f508ddfdbp-1"),
   ("0x1.b274e4p-34", "0x1.9372a6459635fp-2", "0x1.d696173a84a7ap-1"),
   ("0x1.99f79p-34", "0x1.96555b83f5b46p-2", "0x1.d5f717268995bp-1"),
   ("0x1.3a5ep-38", "0x1.993716148d0a5p-2", "0x1.d556f52e7b404p-1"),
   ("0x1.06e368p-34", "0x1.9c17d446c7bcep-2", "0x1.d4b5b1b03af49p-1"),
   ("0x1.79c6ep-36", "0x1.9ef7943cad591p-2", "0x1.d4134d146456ep-1"),
   ("-0x1.014bcp-38", "0x1.a1d6543af46d2p-2", "0x1.d36fc7bcd45bp-1"),
   ("-0x1.38e858p-34", "0x1.a4b41276e9a72p-2", "0x1.d2cb220fa2d8cp-1"),
   ("-0x1.6067cp-34", "0x1.a790cd35df451p-2", "0x1.d2255c70243fep-1"),
   ("-0x1.47f4bp-36", "0x1.aa6c82b4ffa09p-2", "0x1.d17e77444ea0ep-1"),
   ("-0x1.b3e01p-35", "0x1.ad473120f28b9p-2", "0x1.d0d672f6bc322p-1"),
   ("0x1.e2bcp-39", "0x1.b020d6c849efp-2", "0x1.d02d4feb17d94p-1"),
   ("-0x1.75a12p-35", "0x1.b2f971d70af6p-2", "0x1.cf830e8dddb4dp-1"),
   ("-0x1.113d38p-34", "0x1.b5d1009805d13p-2", "0x1.ced7af453b7b4p-1"),
   ("-0x1.5e06p-36", "0x1.b8a7814de51b7p-2", "0x1.ce2b327a1052dp-1"),
   ("0x1.6a6bp-40", "0x1.bb7cf2306be21p-2", "0x1.cd7d9898ab7afp-1"),
   ("-0x1.fcfp-38", "0x1.be51517f482d3p-2", "0x1.cccee20c59781p-1"),
   ("0x1.3ff588p-35", "0x1.c1249d8399431p-2", "0x1.cc1f0f3ef354p-1"),
   ("0x1.23cfdp-35", "0x1.c3f6d47599b2cp-2", "0x1.cb6e209f43591p-1"),
   ("0x1.710cp-37", "0x1.c6c7f49a73b28p-2", "0x1.cabc1699cb32dp-1"),
   ("0x1.1c364p-38", "0x1.c997fc38c9112p-2", "0x1.ca08f19b83549p-1"),
   ("-0x1.ad9d4p-38", "0x1.cc66e9928593fp-2", "0x1.c954b213670c9p-1"),
   ("-0x1.c8cb5p-36", "0x1.cf34baeb9ce5cp-2", "0x1.c89f5870cc0acp-1"),
   ("-0x1.cc923p-36", "0x1.d2016e8c1980bp-2", "0x1.c7e8e522d86d7p-1"),
   ("0x1.06d23p-36", "0x1.d4cd02bbf5125p-2", "0x1.c73158998c2d1p-1"),
   ("-0x1.b56c88p-35", "0x1.d79775b3aa6b7p-2", "0x1.c678b349c3a9dp-1"),
   ("0x1.7b8f5p-36", "0x1.da60c5d1b16dcp-2", "0x1.c5bef59f656bfp-1"),
   ("-0x1.852adp-36", "0x1.dd28f145ffe49p-2", "0x1.c504201344fd9p-1"),
   ("-0x1.11d7p-39", "0x1.dfeff66a649eap-2", "0x1.c44833142899fp-1"),
   ("0x1.238dep-36", "0x1.e2b5d3820348fp-2", "0x1.c38b2f179fe99p-1"),
   ("0x1.5de1ep-36", "0x1.e57a86d5b1691p-2", "0x1.c2cd14929bf71p-1"),
   ("-0x1.9898ap-37", "0x1.e83e0eae6afbap-2", "0x1.c20de3fae39c4p-1"),
   ("-0x1.a22398p-35", "0x1.eb00695aa49eap-2", "0x1.c14d9dc5a0d55p-1"),
   ("0x1.e45e28p-35", "0x1.edc195342ca7fp-2", "0x1.c08c4265b676dp-1"),
   ("0x1.9958dp-35", "0x1.f081907064381p-2", "0x1.bfc9d258e34fp-1"),
   ("-0x1.420a8p-39", "0x1.f3405963c5d13p-2", "0x1.bf064e1546e7fp-1"),
   ("0x1.33ceap-36", "0x1.f5fdee6712454p-2", "0x1.be41b6109ec8cp-1"),
   ("-0x1.5d514p-37", "0x1.f8ba4dbe9af59p-2", "0x1.bd7c0ac73cefdp-1"),
   ("0x1.eb1e2p-37", "0x1.fb7575c39c34fp-2", "0x1.bcb54cb0729f2p-1"),
   ("-0x1.7091p-36", "0x1.fe2f64bc7b28ap-2", "0x1.bbed7c49c847ap-1"),
   ("-0x1.902ap-41", "0x1.00740c82afadep-1", "0x1.bb249a0b712c3p-1"),
   ("-0x1.c359p-35", "0x1.01cfc8725f626p-1", "0x1.ba5aa674be0c3p-1"),
   ("-0x1.be042p-36", "0x1.032ae55dadbdep-1", "0x1.b98fa1fe42a72p-1"),
   ("0x1.26d36p-36", "0x1.0485626ba9785p-1", "0x1.b8c38d26da7cp-1"),
   ("-0x1.56462p-35", "0x1.05df3ec14d8c2p-1", "0x1.b7f6686f8c2b6p-1"),
   ("-0x1.89d9ap-36", "0x1.0738799126acdp-1", "0x1.b728345235eacp-1"),
   ("0x1.ddb76p-36", "0x1.089112046c423p-1", "0x1.b658f14f19e41p-1"),
   ("-0x1.75d84p-36", "0x1.09e9074081746p-1", "0x1.b5889fe9b9be3p-1"),
   ("0x1.a074p-41", "0x1.0b4058790116ap-1", "0x1.b4b7409de23c8p-1"),
   ("0x1.6fb47p-35", "0x1.0c9704d7c454cp-1", "0x1.b3e4d3ee2671cp-1"),
   ("-0x1.67247p-36", "0x1.0ded0b83cc9bap-1", "0x1.b3115a5fcc743p-1"),
   ("0x1.40f08p-38", "0x1.0f426bb2de59ep-1", "0x1.b23cd46fdfd86p-1"),
   ("-0x1.3e38cp-38", "0x1.1097248cd5b15p-1", "0x1.b16742a4eb736p-1"),
   ("-0x1.4907ap-36", "0x1.11eb3540da5e9p-1", "0x1.b090a581da60fp-1"),
   ("0x1.7e9cdp-36", "0x1.133e9cffdfb7fp-1", "0x1.afb8fd8953ef7p-1"),
   ("-0x1.29ea8p-37", "0x1.14915af2d45a5p-1", "0x1.aee04b4400794p-1"),
   ("-0x1.3f98dap-33", "0x1.15e36e4727973p-1", "0x1.ae068f38a0b31p-1"),
   ("-0x1.f08a8p-37", "0x1.1734d63d4a424p-1", "0x1.ad2bc9e287a67p-1"),
   ("0x1.f4d568p-35", "0x1.188591f6368cap-1", "0x1.ac4ffbd240a43p-1"),
   ("0x1.4044fp-35", "0x1.19d5a0a0cf9bcp-1", "0x1.ab7325905720fp-1"),
   ("-0x1.48756p-36", "0x1.1b250170604cfp-1", "0x1.aa9547a35a426p-1"),
   ("-0x1.0b381p-36", "0x1.1c73b39a380bbp-1", "0x1.a9b662915eb3ap-1"),
   ("-0x1.01fdp-39", "0x1.1dc1b64daf83p-1", "0x1.a8d676e553cffp-1"),
   ("-0x1.621c1p-36", "0x1.1f0f08bae2169p-1", "0x1.a7f5852a9a57cp-1"),
   ("0x1.7015d8p-34", "0x1.205baa1b1196bp-1", "0x1.a7138de74aca6p-1"),
   ("0x1.b4708p-36", "0x1.21a79994595c6p-1", "0x1.a63091af6dc26p-1"),
   ("-0x1.b9d398p-34", "0x1.22f2d65e4b185p-1", "0x1.a54c910c24163p-1"),
   ("-0x1.d9b08p-38", "0x1.243d5fb93e635p-1", "0x1.a4678c814ec2fp-1"),
   ("0x1.4f27d8p-35", "0x1.258734cd666bdp-1", "0x1.a38184a465eacp-1"),
   ("-0x1.b2008p-35", "0x1.26d054cba3cefp-1", "0x1.a29a7a05eb03ap-1"),
   ("0x1.73eb6p-37", "0x1.2818bef54af27p-1", "0x1.a1b26d2bb5ff1p-1"),
   ("-0x1.5edp-39", "0x1.296072760dc14p-1", "0x1.a0c95eabc3946p-1"),
   ("0x1.41f19cp-34", "0x1.2aa76e8ae43b8p-1", "0x1.9fdf4f10c6a62p-1"),
   ("-0x1.1ae4d8p-35", "0x1.2bedb25e471a5p-1", "0x1.9ef43ef39f48bp-1"),
   ("-0x1.27d7p-37", "0x1.2d333d348bc93p-1", "0x1.9e082edb869f4p-1"),
   ("0x1.fa6p-38", "0x1.2e780e3ededa4p-1", "0x1.9d1b1f5e6d50ap-1"),
   ("-0x1.1484ep-35", "0x1.2fbc24b2e355fp-1", "0x1.9c2d111009089p-1"),
   ("0x1.51c2ep-35", "0x1.30ff7fcfc1287p-1", "0x1.9b3e047dfc671p-1"),
   ("0x1.ba157p-36", "0x1.32421ec5b0a14p-1", "0x1.9a4dfa41e0baep-1"),
   ("0x1.96c6ep-36", "0x1.338400d1c8557p-1", "0x1.995cf2ecc0ef5p-1"),
   ("0x1.fd47fp-36", "0x1.34c5252d53ef4p-1", "0x1.986aef136621p-1"),
   ("0x1.4f998p-35", "0x1.36058b1209277p-1", "0x1.9777ef4b3e40dp-1"),
   ("0x1.4f0fep-37", "0x1.374531b880713p-1", "0x1.9683f42b87ffdp-1"),
   ("0x1.cde83p-36", "0x1.3884185f1e0e9p-1", "0x1.958efe48096ep-1"),
   ("-0x1.6a026p-34", "0x1.39c23e39e04aep-1", "0x1.94990e3d7d977p-1"),
   ("-0x1.64963p-35", "0x1.3affa2904b79p-1", "0x1.93a2249a7edbp-1"),
   ("0x1.4c462p-36", "0x1.3c3c4498e98ebp-1", "0x1.92aa41fbb951cp-1"),
   ("-0x1.c74bep-36", "0x1.3d78238b3fa7cp-1", "0x1.91b166fe2793ep-1"),
   ("-0x1.45fdcp-37", "0x1.3eb33eab7c36dp-1", "0x1.90b79435c59f5p-1"),
   ("-0x1.1af724p-34", "0x1.3fed95319f619p-1", "0x1.8fbcca4124bacp-1"),
   ("0x1.6b1a9p-36", "0x1.41272664af235p-1", "0x1.8ec109b3d3e33p-1"),
   ("-0x1.25404p-38", "0x1.425ff178b9ff7p-1", "0x1.8dc453318dcdep-1"),
   ("0x1.9f9ap-40", "0x1.4397f5b2b4074p-1", "0x1.8cc6a75177808p-1"),
   ("0x1.acc6ap-35", "0x1.44cf32529a80ap-1", "0x1.8bc806afa62d5p-1"),
   ("0x1.7a45788p-31", "0x1.4605a6af55a26p-1", "0x1.8ac871d63c272p-1"),
   ("0x1.29a3bp-35", "0x1.473b51baeec8dp-1", "0x1.89c7e9a3b27b3p-1"),
   ("-0x1.a889c4p-34", "0x1.4870330209f8cp-1", "0x1.88c66e77d9496p-1"),
   ("0x1.32ecap-36", "0x1.49a449ba6906ap-1", "0x1.87c400fb07b89p-1"),
   ("-0x1.705e1p-36", "0x1.4ad795159561cp-1", "0x1.86c0a1da650c2p-1"),
   ("0x1.40e304p-34", "0x1.4c0a1461bf5dap-1", "0x1.85bc51ac07c91p-1"),
   ("-0x1.1e5018p-35", "0x1.4d3bc6d43485cp-1", "0x1.84b7111c1cf5bp-1"),
   ("-0x1.7028p-38", "0x1.4e6cabbe07a2p-1", "0x1.83b0e0c028ae3p-1"),
   ("0x1.d4a58p-37", "0x1.4f9cc25d55449p-1", "0x1.82a9c13edbbd6p-1"),
   ("-0x1.065fe4p-34", "0x1.50cc09f32d34ap-1", "0x1.81a1b33d75e34p-1"),
   ("-0x1.8054ap-34", "0x1.51fa81ca0eb38p-1", "0x1.8098b75a0237p-1"),
   ("0x1.f320cp-38", "0x1.5328292a7ec46p-1", "0x1.7f8ece35308bcp-1"),
   ("-0x1.bf8aap-35", "0x1.5454ff4f4caa3p-1", "0x1.7e83f87cd6b28p-1"),
   ("0x1.27e92p-36", "0x1.5581038a223bbp-1", "0x1.7d7836cb98d71p-1"),
   ("0x1.4e179cp-34", "0x1.56ac351c8222ep-1", "0x1.7c6b89cb6ebbbp-1"),
   ("-0x1.39219p-36", "0x1.57d693481890ap-1", "0x1.7b5df22750237p-1"),
   ("-0x1.b78c8p-38", "0x1.59001d5f3278ep-1", "0x1.7a4f707c33a4dp-1"),
   ("0x1.8a601p-35", "0x1.5a28d2a79f971p-1", "0x1.79400573528a5p-1"),
   ("0x1.11bcb8p-35", "0x1.5b50b2663331dp-1", "0x1.782fb1b7e7877p-1"),
   ("0x1.9e7978p-35", "0x1.5c77bbe82d185p-1", "0x1.771e75ee7c0a4p-1"),
   ("-0x1.881cep-36", "0x1.5d9dee730248fp-1", "0x1.760c52c3d6c11p-1"),
   ("-0x1.a0d318p-35", "0x1.5ec349565a11ap-1", "0x1.74f948dc4db73p-1"),
   ("-0x1.17fa5p-36", "0x1.5fe7cbddb6e82p-1", "0x1.73e558e110b71p-1"),
   ("0x1.08f7ap-37", "0x1.610b75521e2ap-1", "0x1.72d0837eb839dp-1"),
   ("0x1.f9bp-43", "0x1.622e44fec46d9p-1", "0x1.71bac960e1f67p-1"),
   ("-0x1.e61p-37", "0x1.63503a3138502p-1", "0x1.70a42b31fb4e1p-1"),
   ("-0x1.4f2198p-35", "0x1.647154367b4e7p-1", "0x1.6f8ca99e04332p-1"),
   ("0x1.c83e76p-33", "0x1.659192670b46ap-1", "0x1.6e744546d8881p-1"),
   ("-0x1.91638p-39", "0x1.66b0f3f50f19cp-1", "0x1.6d5afef4c6982p-1"),
   ("-0x1.2f1d8p-35", "0x1.67cf7847c834ap-1", "0x1.6c40d73d66c2p-1"),
   ("0x1.21862p-36", "0x1.68ed1eaabb0f7p-1", "0x1.6b25ced25dddfp-1")]

# global_error(is_sin=true,rel=true)
# i= 2 err= -63.4825073029418?
# global_error(is_sin=false,rel=true)
# i= 216 err= -64.6168214163773?
def global_error(is_sin=true,rel=false):
   global SC
   maxerr = 0
   SC = [(RR(x[0],16),RR(x[1],16),RR(x[2],16)) for x in SC]
   if is_sin:
      for i in range(256):
         errs, sh, sl = evalPSfast ()
         errc, ch, cl = evalPCfast ()
         xi = i/2^11+SC[i][0].exact_rational()
         Si = SC[i][1]
         errSi = abs(n(sin(2*pi*xi)-Si.exact_rational(),200))
         Ci = SC[i][2]
         errCi = abs(n(cos(2*pi*xi)-Ci.exact_rational(),200))
         # s_mul (&sh, &sl, SC[i][2], sh, sl)
         sh_in = sh
         sh, sl = s_mul(Ci, sh, sl)
         err1 = RIFulp(sl)+errCi*sh_in.abs().upper()+Ci*errs
         # s_mul (&ch, &cl, SC[i][1], ch, cl)
         ch_in = ch
         ch, cl = s_mul(Si, ch, cl)
         err2 = RIFulp(cl)+errSi*ch_in.abs().upper()+Si*errc
         # fast_two_sum (h, l, ch, sh)
         h = ch+sh
         u = RIFulp(h)
         l = RIF(-u,u)
         err3 = h.abs().upper()*2^-105
         # *l += sl + cl
         tmp = sl+cl
         l = l+tmp
         err4 = RIFulp(tmp)+RIFulp(l)
         err = err1+err2+err3+err4
         if rel: # convert to relative error
            err =err/(h+l).abs().lower()
         if err>maxerr:
            maxerr = err
            print ("i=", i, "err=", log(err)/log(2.))
   else:
      for i in range(256):
         errs, sh, sl = evalPSfast ()
         errc, ch, cl = evalPCfast ()
         xi = i/2^11+SC[i][0].exact_rational()
         Si = SC[i][1]
         errSi = abs(n(sin(2*pi*xi)-Si.exact_rational(),200))
         Ci = SC[i][2]
         errCi = abs(n(cos(2*pi*xi)-Ci.exact_rational(),200))
         # s_mul (&ch, &cl, SC[i][2], ch, cl)
         ch_in = ch
         ch, cl = s_mul(Ci, ch, cl)
         err1 = RIFulp(cl)+errCi*ch_in.abs().upper()+Ci*errc
         # s_mul (&sh, &sl, SC[i][1], sh, sl
         sh_in = sh
         sh, sl = s_mul(Si, sh, sl)
         err2 = RIFulp(sl)+errSi*sh_in.abs().upper()+Si*errs
         # fast_two_sum (h, l, ch, -sh)
         h = ch-sh
         u = RIFulp(h)
         l = RIF(-u,u)
         err3 = h.abs().upper()*2^-105
         # *l += cl - sl
         tmp = cl - sl
         l = l+tmp
         err4 = RIFulp(tmp)+RIFulp(l)
         err = err1+err2+err3+err4
         if rel: # convert to relative error
            err =err/(h+l).abs().lower()
         if err>maxerr:
            maxerr = err
            print ("i=", i, "err=", log(err)/log(2.))

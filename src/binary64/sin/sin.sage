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
   ("-0x1.6d6f4p-43", "0x1.921f8bea8c466p-9", "0x1.ffff621621d1ep-1"),
   ("0x1.0f66aep-38", "0x1.921f100114f3ep-8", "0x1.fffd8858e8557p-1"),
   ("0x1.9fca9p-41", "0x1.2d96b0e79688bp-7", "0x1.fffa72c978acep-1"),
   ("-0x1.06165cp-40", "0x1.921d1fcab529p-7", "0x1.fff62169b9562p-1"),
   ("-0x1.76633p-39", "0x1.f6a296a2696bp-7", "0x1.fff0943c544d6p-1"),
   ("-0x1.089e2p-42", "0x1.2d936bbdc90a1p-6", "0x1.ffe9cb44b5296p-1"),
   ("0x1.d6fcfap-37", "0x1.5fd4d236c861bp-6", "0x1.ffe1c68708bebp-1"),
   ("-0x1.55cc6p-41", "0x1.92155f792a09fp-6", "0x1.ffd886084d058p-1"),
   ("0x1.dc02fp-39", "0x1.c454f4d42a8fp-6", "0x1.ffce09ce291d4p-1"),
   ("0x1.6c3684p-38", "0x1.f69373260c474p-6", "0x1.ffc251df1b0dfp-1"),
   ("-0x1.0a03cp-38", "0x1.14685db0e8dc1p-5", "0x1.ffb55e42619e1p-1"),
   ("-0x1.8438cep-37", "0x1.2d86574fbfbedp-5", "0x1.ffa72efff5125p-1"),
   ("0x1.67db27p-35", "0x1.46a39722d311bp-5", "0x1.ff97c42075776p-1"),
   ("0x1.982a7cp-38", "0x1.5fc00d2e0dfp-5", "0x1.ff871dadb4adp-1"),
   ("-0x1.6fb75p-40", "0x1.78dbaa5753e8fp-5", "0x1.ff753bb1b9eadp-1"),
   ("0x1.ce4c28p-37", "0x1.91f65f1c32b4fp-5", "0x1.ff621e378def7p-1"),
   ("-0x1.798224p-37", "0x1.ab101bccb776ap-5", "0x1.ff4dc54b23a7ap-1"),
   ("0x1.28bbd2p-36", "0x1.c428d13a98ac8p-5", "0x1.ff3830f8c898bp-1"),
   ("0x1.3ba98ep-35", "0x1.dd406fb6f8e8bp-5", "0x1.ff21614df63bcp-1"),
   ("-0x1.0ccf2p-37", "0x1.f656e798ec418p-5", "0x1.ff095658ed93fp-1"),
   ("0x1.87fddp-38", "0x1.07b614e6c97c5p-4", "0x1.fef0102821248p-1"),
   ("-0x1.04b098p-38", "0x1.1440134bd80c4p-4", "0x1.fed58ecb6abp-1"),
   ("-0x1.4d956dp-35", "0x1.20c9673e7ec3ap-4", "0x1.feb9d25329029p-1"),
   ("-0x1.33798cp-36", "0x1.2d5209255aefap-4", "0x1.fe9cdad02a478p-1"),
   ("-0x1.f1dep-38", "0x1.39d9f1294e69ep-4", "0x1.fe7ea8548a539p-1"),
   ("0x1.bb85ap-39", "0x1.46611793cd453p-4", "0x1.fe5f3af2e01bdp-1"),
   ("-0x1.4db84cp-37", "0x1.52e774a0bfff7p-4", "0x1.fe3e92bea8604p-1"),
   ("0x1.c4d0a6p-36", "0x1.5f6d00b4bce2cp-4", "0x1.fe1cafcbb72d4p-1"),
   ("-0x1.5d3afdp-35", "0x1.6bf1b3d687df2p-4", "0x1.fdf9922fa3ee8p-1"),
   ("0x1.8c51dcp-33", "0x1.787586f35294ep-4", "0x1.fdd539fe3a676p-1"),
   ("-0x1.09a3acp-37", "0x1.84f87128d447bp-4", "0x1.fdafa7514f20cp-1"),
   ("0x1.db96c2p-35", "0x1.917a6bd9d6e8bp-4", "0x1.fd88da3cc919p-1"),
   ("0x1.dc86a4p-36", "0x1.9dfb6ebded4aap-4", "0x1.fd60d2da4ff64p-1"),
   ("-0x1.4c3f9p-37", "0x1.aa7b724087a2bp-4", "0x1.fd3791422e4e1p-1"),
   ("0x1.ab769ep-35", "0x1.b6fa6ed86c1e4p-4", "0x1.fd0d158d3e124p-1"),
   ("-0x1.014fp-41", "0x1.c3785c79b9f66p-4", "0x1.fce15fd6db19ep-1"),
   ("-0x1.1ede42p-35", "0x1.cff533a50a2b1p-4", "0x1.fcb47039473fdp-1"),
   ("0x1.1406fcp-36", "0x1.dc70ecc1a48b6p-4", "0x1.fc8646cfd23a9p-1"),
   ("0x1.1e7decp-34", "0x1.e8eb7ffa377ccp-4", "0x1.fc56e3b76e406p-1"),
   ("-0x1.711e6p-36", "0x1.f564e56199645p-4", "0x1.fc26470e3d7ap-1"),
   ("-0x1.0548ap-39", "0x1.00ee8ad695ba2p-3", "0x1.fbf470f0ac105p-1"),
   ("0x1.42685cp-35", "0x1.072a04838125dp-3", "0x1.fbc1617e0304ep-1"),
   ("0x1.40cc18p-37", "0x1.0d64dbcd19ffbp-3", "0x1.fb8d18d65a494p-1"),
   ("-0x1.acb79p-36", "0x1.139f0ce878bd7p-3", "0x1.fb5797198ac4cp-1"),
   ("-0x1.e70ccp-39", "0x1.19d8940b24dccp-3", "0x1.fb20dc6823e96p-1"),
   ("-0x1.67e3bp-38", "0x1.20116d4dafe44p-3", "0x1.fae8e8e476ec3p-1"),
   ("0x1.9ba8p-36", "0x1.264994e4d3165p-3", "0x1.faafbcb0a1685p-1"),
   ("-0x1.e4b828p-36", "0x1.2c8106e303c49p-3", "0x1.fa7557f0c22d9p-1"),
   ("0x1.0f233ep-34", "0x1.32b7bfa17a346p-3", "0x1.fa39bac721e72p-1"),
   ("-0x1.e9fb32p-35", "0x1.38edbb00f6534p-3", "0x1.f9fce55b50c6ap-1"),
   ("-0x1.230a2p-39", "0x1.3f22f57d43a42p-3", "0x1.f9bed7cfc2566p-1"),
   ("0x1.96564p-36", "0x1.45576b17803b5p-3", "0x1.f97f924c5de72p-1"),
   ("0x1.03f8b2p-35", "0x1.4b8b17fdeb8fap-3", "0x1.f93f14f818a4cp-1"),
   ("0x1.dff6f8p-36", "0x1.51bdf85f4b947p-3", "0x1.f8fd5ffaa5f32p-1"),
   ("-0x1.f20a6cp-36", "0x1.57f0085f46529p-3", "0x1.f8ba737cf6685p-1"),
   ("-0x1.61388p-41", "0x1.5e21444891d1bp-3", "0x1.f8764fa71636p-1"),
   ("0x1.36a566p-34", "0x1.6451a840dc58fp-3", "0x1.f830f4a362954p-1"),
   ("0x1.3e5bf8p-36", "0x1.6a8130533d08p-3", "0x1.f7ea629e3794ep-1"),
   ("0x1.2f79ep-38", "0x1.70afd8d176c49p-3", "0x1.f7a299c198688p-1"),
   ("0x1.ef77cp-40", "0x1.76dd9de56b974p-3", "0x1.f7599a3a0d93dp-1"),
   ("0x1.ada62p-39", "0x1.7d0a7bbdd2789p-3", "0x1.f70f6434b0126p-1"),
   ("0x1.3ae8ep-36", "0x1.83366e8d91c5dp-3", "0x1.f6c3f7df2cf85p-1"),
   ("0x1.cff1d8p-36", "0x1.89617281d7e2dp-3", "0x1.f67755683dd15p-1"),
   ("0x1.5cc3dp-36", "0x1.8f8b83cacd007p-3", "0x1.f6297cff405aap-1"),
   ("-0x1.6eac5cp-34", "0x1.95b49e89be3fep-3", "0x1.f5da6ed51ab84p-1"),
   ("0x1.920a78p-35", "0x1.9bdcbf376eb38p-3", "0x1.f58a2b170ae7ap-1"),
   ("-0x1.7ea3ap-34", "0x1.a203e19f1ff4ap-3", "0x1.f538b1fbe82c4p-1"),
   ("0x1.d656cp-38", "0x1.a82a025c69a9dp-3", "0x1.f4e603b09fd25p-1"),
   ("0x1.9ef77p-36", "0x1.ae4f1d643628cp-3", "0x1.f492206b8640dp-1"),
   ("0x1.5f1ad4p-35", "0x1.b4732efc41d34p-3", "0x1.f43d085f83a54p-1"),
   ("0x1.c4fc7p-36", "0x1.ba9633548353p-3", "0x1.f3e6bbc16ee3bp-1"),
   ("0x1.65ef6cp-33", "0x1.c0b826ca2e464p-3", "0x1.f38f3ac46197dp-1"),
   ("-0x1.10d4e3p-33", "0x1.c6d9051bb98cbp-3", "0x1.f33685a527a87p-1"),
   ("0x1.3fc51p-36", "0x1.ccf8cb34fdf74p-3", "0x1.f2dc9c9051226p-1"),
   ("-0x1.2f7aap-33", "0x1.d31774b5c9a23p-3", "0x1.f2817fc61381p-1"),
   ("-0x1.931ba3p-33", "0x1.d934fe2dd3228p-3", "0x1.f2252f79ace71p-1"),
   ("0x1.69c1e6p-34", "0x1.df5164015451fp-3", "0x1.f1c7abe17a73cp-1"),
   ("-0x1.02197dp-33", "0x1.e56ca1c863babp-3", "0x1.f168f540f266ap-1"),
   ("-0x1.2c34ep-38", "0x1.eb86b461f9509p-3", "0x1.f1090bc8a71b9p-1"),
   ("0x1.2e9ad9p-33", "0x1.f19f97cee7502p-3", "0x1.f0a7efb75511ap-1"),
   ("0x1.6977fp-36", "0x1.f7b7481020341p-3", "0x1.f045a14cb1655p-1"),
   ("0x1.26f354p-35", "0x1.fdcdc1b501c42p-3", "0x1.efe220c0460a5p-1"),
   ("-0x1.7bcc98p-36", "0x1.01f180695e83ap-2", "0x1.ef7d6e52155fbp-1"),
   ("-0x1.4dba3cp-35", "0x1.04fb80df8a0adp-2", "0x1.ef178a3eccd71p-1"),
   ("-0x1.037653p-32", "0x1.0804e04619dc6p-2", "0x1.eeb074c852fcbp-1"),
   ("-0x1.32e8d8p-35", "0x1.0b0d9cfa1ae77p-2", "0x1.ee482e262795dp-1"),
   ("0x1.042ccp-37", "0x1.0e15b4e239b7cp-2", "0x1.eddeb6a05d726p-1"),
   ("0x1.84f19p-37", "0x1.111d262c45d07p-2", "0x1.ed740e765bd98p-1"),
   ("0x1.24ce59p-33", "0x1.1423ef0a40699p-2", "0x1.ed0835e7a8e0fp-1"),
   ("-0x1.d144p-38", "0x1.172a0d76b54d9p-2", "0x1.ec9b2d3c54ep-1"),
   ("0x1.529b1p-37", "0x1.1a2f7fbf8ec8ep-2", "0x1.ec2cf4b18ac69p-1"),
   ("0x1.9d211p-33", "0x1.1d34440847cf7p-2", "0x1.ebbd8c8b1dbe6p-1"),
   ("0x1.7df932p-34", "0x1.203858467177p-2", "0x1.eb4cf51466bebp-1"),
   ("-0x1.643948p-34", "0x1.233bbab3d9ed1p-2", "0x1.eadb2e8fb8d16p-1"),
   ("-0x1.3800dcp-35", "0x1.263e6991aa723p-2", "0x1.ea68393ef253cp-1"),
   ("0x1.86d7f4p-34", "0x1.294062f687f02p-2", "0x1.e9f4156afe6fp-1"),
   ("-0x1.46a5acp-35", "0x1.2c41a4e57f3bfp-2", "0x1.e97ec360ad262p-1"),
   ("0x1.30db5p-36", "0x1.2f422db0899b5p-2", "0x1.e908436198969p-1"),
   ("0x1.5a13ccp-35", "0x1.3241fb6799231p-2", "0x1.e89095ba336c9p-1"),
   ("-0x1.a8d808p-36", "0x1.35410c2b9be67p-2", "0x1.e817bab531d62p-1"),
   ("-0x1.16b238p-33", "0x1.383f5e2834027p-2", "0x1.e79db29c675bcp-1"),
   ("0x1.037ef4p-35", "0x1.3b3cefa348efp-2", "0x1.e7227db62bf85p-1"),
   ("0x1.8ebb1p-35", "0x1.3e39be9b92c75p-2", "0x1.e6a61c5512965p-1"),
   ("0x1.7775p-36", "0x1.4135c943a6601p-2", "0x1.e6288ec431914p-1"),
   ("-0x1.68b468p-36", "0x1.44310dc679fc2p-2", "0x1.e5a9d550a02ddp-1"),
   ("-0x1.13721cp-35", "0x1.472b8a523d09dp-2", "0x1.e529f047b43cep-1"),
   ("0x1.91511p-37", "0x1.4a253d12e28bfp-2", "0x1.e4a8dff7ea167p-1"),
   ("0x1.481e24p-34", "0x1.4d1e242f2bf36p-2", "0x1.e426a4b16cc27p-1"),
   ("0x1.c007dp-37", "0x1.50163dc2e3684p-2", "0x1.e3a33ec723297p-1"),
   ("-0x1.12634p-37", "0x1.530d880a28693p-2", "0x1.e31eae87308fbp-1"),
   ("0x1.552092p-34", "0x1.560401372ac2bp-2", "0x1.e298f4422ba69p-1"),
   ("-0x1.2c7516p-33", "0x1.58f9a74ccf0aep-2", "0x1.e2121051e46bep-1"),
   ("-0x1.c3428ep-33", "0x1.5bee78a505e77p-2", "0x1.e18a030189d08p-1"),
   ("0x1.30e4e8p-36", "0x1.5ee2737bac57dp-2", "0x1.e100cca245fccp-1"),
   ("-0x1.a42e4ap-34", "0x1.61d595bedeabcp-2", "0x1.e0766d944915ep-1"),
   ("0x1.84d1a2p-34", "0x1.64c7dddce46ccp-2", "0x1.dfeae6213249dp-1"),
   ("0x1.1e77c8p-34", "0x1.67b949d16b735p-2", "0x1.df5e36a87e33p-1"),
   ("0x1.b12cp-43", "0x1.6aa9d7dc7cda1p-2", "0x1.ded05f7de38cap-1"),
   ("0x1.7e4744p-35", "0x1.6d99863d029aep-2", "0x1.de4160f602742p-1"),
   ("-0x1.b3f248p-36", "0x1.7088530d25742p-2", "0x1.ddb13b6d475d6p-1"),
   ("-0x1.6adecep-34", "0x1.73763c8a145b5p-2", "0x1.dd1fef3a469f6p-1"),
   ("-0x1.89a718p-34", "0x1.766340e943686p-2", "0x1.dc8d7cb5d44dp-1"),
   ("0x1.13fdebp-33", "0x1.794f5e6dd62cfp-2", "0x1.dbf9e436d8639p-1"),
   ("0x1.2e56b8p-35", "0x1.7c3a93154eb99p-2", "0x1.db652622d9b1dp-1"),
   ("0x1.a7e4a4p-34", "0x1.7f24dd40da0ep-2", "0x1.dacf42cc76656p-1"),
   ("-0x1.163d6p-38", "0x1.820e3b0485789p-2", "0x1.da383a967d314p-1"),
   ("0x1.cae3e4p-34", "0x1.84f6aab9a433ep-2", "0x1.d9a00dd690397p-1"),
   ("-0x1.30f98p-35", "0x1.87de2a6775698p-2", "0x1.d906bcf3e027cp-1"),
   ("0x1.7fcb88p-35", "0x1.8ac4b871b75bp-2", "0x1.d86c484371db1p-1"),
   ("-0x1.3e10ccp-34", "0x1.8daa52e558adep-2", "0x1.d7d0b02d12dbap-1"),
   ("0x1.c98e98p-35", "0x1.908ef82422a84p-2", "0x1.d733f507a7babp-1"),
   ("0x1.077686p-32", "0x1.9372a6538f321p-2", "0x1.d696173785fd9p-1"),
   ("0x1.3568bfp-32", "0x1.96555b969b284p-2", "0x1.d5f7172281cc9p-1"),
   ("0x1.b01174p-34", "0x1.9937161dd4714p-2", "0x1.d556f52c7579bp-1"),
   ("-0x1.68a6dp-33", "0x1.9c17d430aabbbp-2", "0x1.d4b5b1b517413p-1"),
   ("0x1.03c574p-33", "0x1.9ef79446372d3p-2", "0x1.d4134d1247253p-1"),
   ("0x1.26ec78p-35", "0x1.a1d6543e9e8f1p-2", "0x1.d36fc7bc02b56p-1"),
   ("-0x1.ed9d1ep-33", "0x1.a4b41267d2ddep-2", "0x1.d2cb221309315p-1"),
   ("-0x1.4e9a408p-31", "0x1.a790cd01ee851p-2", "0x1.d2255c7bf0c9ap-1"),
   ("-0x1.19827cp-32", "0x1.aa6c829db385dp-2", "0x1.d17e7749a4824p-1"),
   ("-0x1.b3e01p-35", "0x1.ad473120f28b9p-2", "0x1.d0d672f6bc322p-1"),
   ("0x1.11144p-36", "0x1.b020d6c978e35p-2", "0x1.d02d4fead156ep-1"),
   ("-0x1.75a12p-35", "0x1.b2f971d70af6p-2", "0x1.cf830e8dddb4dp-1"),
   ("-0x1.99e604p-33", "0x1.b5d1008be567dp-2", "0x1.ced7af48199ecp-1"),
   ("-0x1.5e06p-36", "0x1.b8a7814de51b7p-2", "0x1.ce2b327a1052dp-1"),
   ("0x1.c6e44p-37", "0x1.bb7cf2318dd6bp-2", "0x1.cd7d989865d1bp-1"),
   ("0x1.fad144p-33", "0x1.be515196601fep-2", "0x1.cccee206c1f56p-1"),
   ("0x1.3ff588p-35", "0x1.c1249d8399431p-2", "0x1.cc1f0f3ef354p-1"),
   ("-0x1.2c1eeap-33", "0x1.c3f6d4652ae08p-2", "0x1.cb6e20a34df66p-1"),
   ("-0x1.5e52bp-34", "0x1.c6c7f491bbdap-2", "0x1.cabc169bf459ep-1"),
   ("0x1.1c364p-38", "0x1.c997fc38c9112p-2", "0x1.ca08f19b83549p-1"),
   ("0x1.694adp-35", "0x1.cc66e997121cbp-2", "0x1.c954b21241f5ep-1"),
   ("-0x1.c8cb5p-36", "0x1.cf34baeb9ce5cp-2", "0x1.c89f5870cc0acp-1"),
   ("0x1.fdd48p-36", "0x1.d2016e9166d0ep-2", "0x1.c7e8e5217d95cp-1"),
   ("0x1.3f9ab4p-34", "0x1.d4cd02c17f5dap-2", "0x1.c73158981f00ap-1"),
   ("0x1.0102e4p-34", "0x1.d79775be07a02p-2", "0x1.c678b347136p-1"),
   ("-0x1.161284p-33", "0x1.da60c5c38845cp-2", "0x1.c5bef5a318eb9p-1"),
   ("0x1.0e8a0cp-34", "0x1.dd28f14dfccadp-2", "0x1.c50420112a85p-1"),
   ("-0x1.8bbap-37", "0x1.dfeff66981909p-2", "0x1.c448331464d6p-1"),
   ("0x1.4db844p-32", "0x1.e2b5d39d54525p-2", "0x1.c38b2f1052fcp-1"),
   ("-0x1.de571p-36", "0x1.e57a86d137f2p-2", "0x1.c2cd1493d05c2p-1"),
   ("-0x1.9898ap-37", "0x1.e83e0eae6afbap-2", "0x1.c20de3fae39c4p-1"),
   ("-0x1.a22398p-35", "0x1.eb00695aa49eap-2", "0x1.c14d9dc5a0d55p-1"),
   ("0x1.e45e28p-35", "0x1.edc195342ca7fp-2", "0x1.c08c4265b676dp-1"),
   ("0x1.dfd6ep-34", "0x1.f08190764c4edp-2", "0x1.bfc9d2574028bp-1"),
   ("0x1.c77014p-33", "0x1.f340597781ec3p-2", "0x1.bf064e0fc4516p-1"),
   ("0x1.33ceap-36", "0x1.f5fdee6712454p-2", "0x1.be41b6109ec8cp-1"),
   ("-0x1.3f3bc4p-32", "0x1.f8ba4da444c73p-2", "0x1.bd7c0aceb2a2dp-1"),
   ("-0x1.35b148p-34", "0x1.fb7575bbb310fp-2", "0x1.bcb54cb2b458ap-1"),
   ("-0x1.7091p-36", "0x1.fe2f64bc7b28ap-2", "0x1.bbed7c49c847ap-1"),
   ("-0x1.f6c96p-36", "0x1.00740c8162663p-1", "0x1.bb249a0c320bfp-1"),
   ("-0x1.370108p-33", "0x1.01cfc86e2ba24p-1", "0x1.ba5aa67731038p-1"),
   ("0x1.b9bcbp-35", "0x1.032ae56132445p-1", "0x1.b98fa1fc321a7p-1"),
   ("0x1.2b7fap-36", "0x1.0485626baca12p-1", "0x1.b8c38d26d89dfp-1"),
   ("-0x1.56462p-35", "0x1.05df3ec14d8c2p-1", "0x1.b7f6686f8c2b6p-1"),
   ("0x1.ba509p-35", "0x1.0738799483efp-1", "0x1.b728345031b2dp-1"),
   ("-0x1.6a3b6cp-34", "0x1.089111ff5cc03p-1", "0x1.b658f15227cdp-1"),
   ("0x1.3d347ep-33", "0x1.09e90748238d1p-1", "0x1.b5889fe51625p-1"),
   ("-0x1.1611ap-32", "0x1.0b40586d53db2p-1", "0x1.b4b740a50784bp-1"),
   ("0x1.5e1834p-34", "0x1.0c9704d980f6dp-1", "0x1.b3e4d3ed14783p-1"),
   ("0x1.baa1acp-32", "0x1.0ded0b9732cap-1", "0x1.b3115a53c3527p-1"),
   ("0x1.40f08p-38", "0x1.0f426bb2de59ep-1", "0x1.b23cd46fdfd86p-1"),
   ("0x1.0efdap-36", "0x1.1097248dbebf4p-1", "0x1.b16742a458dedp-1"),
   ("0x1.f1ae6p-35", "0x1.11eb3544492a4p-1", "0x1.b090a57fade33p-1"),
   ("-0x1.1e467ep-33", "0x1.133e9cf8f5a2ap-1", "0x1.afb8fd8dbc738p-1"),
   ("0x1.bf509p-36", "0x1.14915af45e767p-1", "0x1.aee04b4303815p-1"),
   ("-0x1.f6faacp-33", "0x1.15e36e435fd85p-1", "0x1.ae068f3b1211ep-1"),
   ("-0x1.68fefcp-34", "0x1.1734d63a37129p-1", "0x1.ad2bc9e487c56p-1"),
   ("-0x1.e81af4p-32", "0x1.188591df9949cp-1", "0x1.ac4ffbe1104c5p-1"),
   ("0x1.ae5904p-34", "0x1.19d5a0a39452bp-1", "0x1.ab73258e83d86p-1"),
   ("-0x1.1c7d84p-33", "0x1.1b25016b65f02p-1", "0x1.aa9547a6a81e1p-1"),
   ("-0x1.bd9798p-34", "0x1.1c73b3965a99bp-1", "0x1.a9b66293f3d3fp-1"),
   ("0x1.0967a08p-31", "0x1.1dc1b66363536p-1", "0x1.a8d676d6bad49p-1"),
   ("0x1.f95e8p-36", "0x1.1f0f08bd110bcp-1", "0x1.a7f585291fe0bp-1"),
   ("-0x1.392476p-33", "0x1.205baa10fc3d4p-1", "0x1.a7138dee2a2c2p-1"),
   ("-0x1.0fd451p-32", "0x1.21a799883dfcdp-1", "0x1.a63091b7bc2b1p-1"),
   ("-0x1.b9d398p-34", "0x1.22f2d65e4b185p-1", "0x1.a54c910c24163p-1"),
   ("-0x1.d9b08p-38", "0x1.243d5fb93e635p-1", "0x1.a4678c814ec2fp-1"),
   ("-0x1.f1b45cp-34", "0x1.258734c6b5f27p-1", "0x1.a38184a914227p-1"),
   ("-0x1.997874p-34", "0x1.26d054c9b5726p-1", "0x1.a29a7a07472edp-1"),
   ("0x1.73eb6p-37", "0x1.2818bef54af27p-1", "0x1.a1b26d2bb5ff1p-1"),
   ("-0x1.5edp-39", "0x1.296072760dc14p-1", "0x1.a0c95eabc3946p-1"),
   ("0x1.666cecp-33", "0x1.2aa76e8ed3f0ep-1", "0x1.9fdf4f0df2f97p-1"),
   ("0x1.0e955ep-32", "0x1.2bedb26a7300ep-1", "0x1.9ef43eead31dp-1"),
   ("0x1.55246p-37", "0x1.2d333d355610ap-1", "0x1.9e082edaf377dp-1"),
   ("0x1.eaf96p-37", "0x1.2e780e3f2a31ep-1", "0x1.9d1b1f5e36269p-1"),
   ("-0x1.1484ep-35", "0x1.2fbc24b2e355fp-1", "0x1.9c2d111009089p-1"),
   ("0x1.6aa5c4p-34", "0x1.30ff7fd1aa199p-1", "0x1.9b3e047c91c75p-1"),
   ("0x1.ba157p-36", "0x1.32421ec5b0a14p-1", "0x1.9a4dfa41e0baep-1"),
   ("0x1.96c6ep-36", "0x1.338400d1c8557p-1", "0x1.995cf2ecc0ef5p-1"),
   ("0x1.4303ccp-33", "0x1.34c5253267d4ap-1", "0x1.986aef0f8f633p-1"),
   ("-0x1.6289cp-34", "0x1.36058b0cef351p-1", "0x1.9777ef4f1fe16p-1"),
   ("0x1.058624p-34", "0x1.374531baa44dcp-1", "0x1.9683f429e47ffp-1"),
   ("-0x1.59e34p-35", "0x1.3884185c50544p-1", "0x1.958efe4a327f6p-1"),
   ("0x1.90643p-34", "0x1.39c23e414503ap-1", "0x1.94990e37c1d17p-1"),
   ("-0x1.099ef4p-34", "0x1.3affa28f73311p-1", "0x1.93a2249b27a49p-1"),
   ("0x1.4c462p-36", "0x1.3c3c4498e98ebp-1", "0x1.92aa41fbb951cp-1"),
   ("-0x1.c74bep-36", "0x1.3d78238b3fa7cp-1", "0x1.91b166fe2793ep-1"),
   ("-0x1.5299e8p-34", "0x1.3eb33ea89fdd3p-1", "0x1.90b794380c141p-1"),
   ("-0x1.70ff52p-33", "0x1.3fed952d434dep-1", "0x1.8fbcca44a1fap-1"),
   ("0x1.6b1a9p-36", "0x1.41272664af235p-1", "0x1.8ec109b3d3e33p-1"),
   ("-0x1.25404p-38", "0x1.425ff178b9ff7p-1", "0x1.8dc453318dcdep-1"),
   ("0x1.471529p-32", "0x1.4397f5bf1576ap-1", "0x1.8cc6a7475ea55p-1"),
   ("0x1.acc6ap-35", "0x1.44cf32529a80ap-1", "0x1.8bc806afa62d5p-1"),
   ("0x1.8f051ap-31", "0x1.4605a6b0e7b7cp-1", "0x1.8ac871d4f01a1p-1"),
   ("0x1.29a3bp-35", "0x1.473b51baeec8dp-1", "0x1.89c7e9a3b27b3p-1"),
   ("-0x1.bcdb44p-34", "0x1.48703301d900fp-1", "0x1.88c66e78023bdp-1"),
   ("-0x1.137d28p-35", "0x1.49a449b86575dp-1", "0x1.87c400fcb988p-1"),
   ("-0x1.705e1p-36", "0x1.4ad795159561cp-1", "0x1.86c0a1da650c2p-1"),
   ("0x1.40e304p-34", "0x1.4c0a1461bf5dap-1", "0x1.85bc51ac07c91p-1"),
   ("-0x1.543d64p-33", "0x1.4d3bc6cf32f06p-1", "0x1.84b71120679f3p-1"),
   ("-0x1.7028p-38", "0x1.4e6cabbe07a2p-1", "0x1.83b0e0c028ae3p-1"),
   ("0x1.d4a58p-37", "0x1.4f9cc25d55449p-1", "0x1.82a9c13edbbd6p-1"),
   ("0x1.1ef05p-33", "0x1.50cc09fae7f34p-1", "0x1.81a1b336b5b8ep-1"),
   ("0x1.5c62b2p-33", "0x1.51fa81d405f1bp-1", "0x1.8098b75140371p-1"),
   ("-0x1.9f512cp-33", "0x1.5328292292779p-1", "0x1.7f8ece3c320c6p-1"),
   ("-0x1.bdbec8p-32", "0x1.5454ff4101109p-1", "0x1.7e83f8898eb84p-1"),
   ("0x1.a5ef5ap-33", "0x1.558103912c4bap-1", "0x1.7d7836c54b935p-1"),
   ("0x1.038fa7p-32", "0x1.56ac3522edc9ep-1", "0x1.7c6b89c5a62c5p-1"),
   ("-0x1.267bbcp-33", "0x1.57d6934373cfep-1", "0x1.7b5df22b858fbp-1"),
   ("0x1.44f8p-34", "0x1.59001d6264962p-1", "0x1.7a4f707949905p-1"),
   ("0x1.992d04p-34", "0x1.5a28d2a98a4bbp-1", "0x1.7940057190469p-1"),
   ("0x1.11bcb8p-35", "0x1.5b50b2663331dp-1", "0x1.782fb1b7e7877p-1"),
   ("0x1.9e7978p-35", "0x1.5c77bbe82d185p-1", "0x1.771e75ee7c0a4p-1"),
   ("0x1.1b40dp-32", "0x1.5d9dee7e0baf4p-1", "0x1.760c52b985e5cp-1"),
   ("-0x1.a0d318p-35", "0x1.5ec349565a11ap-1", "0x1.74f948dc4db73p-1"),
   ("-0x1.41c65p-34", "0x1.5fe7cbdb785dp-1", "0x1.73e558e330601p-1"),
   ("0x1.48886p-36", "0x1.610b75528dae6p-1", "0x1.72d0837e4e0d8p-1"),
   ("-0x1.2e2d9p-34", "0x1.622e44fc14a81p-1", "0x1.71bac96374cf6p-1"),
   ("-0x1.e61p-37", "0x1.63503a3138502p-1", "0x1.70a42b31fb4e1p-1"),
   ("-0x1.9793p-35", "0x1.64715436299e2p-1", "0x1.6f8ca99e536b8p-1"),
   ("-0x1.df04bep-33", "0x1.659192569d55ep-1", "0x1.6e744556e07d3p-1"),
   ("-0x1.0fc07p-35", "0x1.66b0f3f3fa9d6p-1", "0x1.6d5afef5d6097p-1"),
   ("0x1.0efc14p-33", "0x1.67cf784dd6426p-1", "0x1.6c40d7376b9c4p-1"),
   ("0x1.1ce5e4p-34", "0x1.68ed1eac9499cp-1", "0x1.6b25ced087393p-1"),
]

# global_error(is_sin=true,rel=true)
# i= 2 err= -65.0210177185380?
# global_error(is_sin=false,rel=true)
# i= 253 err= -65.5249308469996?
def global_error(is_sin=true,rel=false):
   global SC
   maxerr = 0
   SC = [(RR(x[0],16),RR(x[1],16),RR(x[2],16)) for x in SC]
   # check SC[i][0] is integer multiple of 2^-62, and |SC[i][0]| < 2^-30
   for i in range(256):
      x = SC[i][0]
      X = x.exact_rational()
      assert X/2^-62 in ZZ, "X/2^-62 in ZZ"
      assert abs(X)<2^-30, "abs(X)<2^-30"
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

from functools import cmp_to_key

# entries are (t0,t1,e)
def cmp(x,y):
   xmin = x[0]*2^x[2]
   xmax = x[1]*2^x[2]
   ymin = y[0]*2^y[2]
   ymax = y[1]*2^y[2]
   if xmax <= ymin:
      return int(-1)
   if ymax <= xmin:
      return int(1)
   if (xmin <= ymin and xmax < ymax):
      return int(-1)
   if (xmin < ymin and xmax <= ymax):
      return int(-1)
   if (ymin <= xmin and ymax < xmax):
      return int(1)
   if (ymin < xmin and ymax <= xmax):
      return int(1)
   return int(0)

# statall("out")
# [((1, -1074), (4503599627370496, -1074)),
# [((4503599627370496, -113), (6667114472105511, -50))]
def statall(f):
   f = open(f,"r")
   l = []
   while true:
      s = f.readline()
      if s=='':
         break
      s = s.split(" ")
      assert len(s) == 5
      t0 = ZZ(s[0])
      t1 = ZZ(s[1])
      e = ZZ(s[2])
      n = ZZ(s[3])
      nn = ZZ(s[4])
      assert nn == 53, "nn == 53"
      assert t0.nbits() == n, "t0.nbits() == n"
      assert (t1-1).nbits() == n, "(t1-1).nbits() == n"
      l.append((t0,t1,e-n))
   f.close()
   l.sort(key=cmp_to_key(cmp))
   l2 = []
   for t0,t1,e in l:
      if l2==[]:
         l2 = [((t0,e),(t1,e))]
      else:
         t1old,e1old = l2[-1][1]
         if t1old*2^e1old > t0*2^e:
            print ((t1old,e1old), (t0, e))
         assert t1old*2^e1old <= t0*2^e, "t1old*2^e1old <= t0*2^e"
         if t1old*2^e1old == t0*2^e:
            l2[-1] = (l2[-1][0],(t1,e))
         else:
            l2.append(((t0,e),(t1,e)))
   l = l2
   return l

# return the number of terms needed to evaluate erf(x) with prec bits of
# accuracy
# erf x = 2/sqrt(pi) sum((-1)^n*x^(2n+1)/n!/(2n+1),n=0..infinity)
# nb_terms_taylor(RR(1),70)
# 20
# nb_terms_taylor(RR(1/2),70)
# 14
# nb_terms_taylor(RR(1/4),70)
# 10
# nb_terms_taylor(5,70)
# 102
def nb_terms_taylor(x,prec):
   x = RR(x)
   y = erf(x)
   eps = abs(y)/2^prec
   t = numerical_approx(2/sqrt(pi)*x)
   n = 0
   while true:
      n += 1
      t = t*x^2/n
      if abs(t/(2*n+1))<eps:
         return n-1

# same for the asymptotic expansion
# erfc(x) = exp(-x^2)/x/sqrt(pi)*(1+sum((-1)^n*(1*3*5*...*(2n-1))/(2x^2)^n,n=1..infinity))
# nb_terms_asympt(5,70)
# 13
# nb_terms_asympt(4.9,70)
# 16
def nb_terms_asympt(x,prec):
   x = RR(x)
   y = abs(erf(x))
   eps = y/2^prec
   t = mint = abs(numerical_approx(exp(-x^2)/x/sqrt(pi)))
   n = 0
   while true:
      n += 1
      oldt = t
      t = abs(t*(2*n-1)/(2*x^2))
      if t < mint:
         mint = t
      if t<eps:
         return n-1
      if t>oldt:
         return "does not converge", mint, log(y/mint)

def double_fact(n):
   return product(k for k in range(1,n+1,2))

# return a minimax polynomial for erf(x0+x) over -h <= x <= h, of degree d
# with coefficients of degree 0..i of precision 107,
# and larger coefficients of precision 53
# sollya_approx(RR(0.2),RR(0.1),6,2)
# (-0.03403920444096596165106660691890*x^6 + 0.09095324353121236882202538254205*x^5 + 0.1055211058400329360296154845855*x^4 - 0.3324662762010896499553780358838*x^3 - 0.2168269541858736853519732891500*x^2 + 1.084134784949880605455645490444*x + 0.2227025892091853463730483303844, 34.9180000000000)
# sollya_approx(RR(0.2),RR(0.1),6,2,rel=true)
# (-0.03352012072306500417573360550705*x^6 + 0.09095692739678004623549156804074*x^5 + 0.1055133196294501252143760439139*x^4 - 0.3324663222575241494638476069667*x^3 - 0.2168269249878638063019510679249*x^2 + 1.084134785065062106344918858833*x + 0.2227025891929645116318119446500, 32.6330000000000)
def sollya_approx(x0,h,d,i,rel=false):
   if rel==false:
      s = "./doit.sh"
   else:
      s = "./doit_rel.sh"
   s += " " + get_hex(x0) + " " + get_hex(h) + " " + str(d)
   s += " '107"
   for j in range(i):
      s += ",107"
   s += ",53...'"
   os.system (s)
   f = open("out","r")
   s = f.readline()
   R.<x> = RealField(107)[]
   p = R(s)
   s = f.readline()
   ret = RR(s[:-1])
   f.close()
   return p, ret

# return the worst approximation for erf(j0*2*h+x), -h <= x <= h, j0 <= j < j1
# sollya_approx_all(1,63,1/32,10,1)
# j= 14 err= 70.2210000000000
# sollya_approx_all(1,95,1/32,10,1,rel=true)
# j= 1 err= 68.9360000000000
def sollya_approx_all(j0,j1,h,d,i,rel=false):
   wo = infinity
   h = RR(h)
   for j in range(j0,j1):
      _, err = sollya_approx(j*2*h,h,d,i,rel=rel)
      if err<wo:
         print ("j=", j, "err=", err)
         wo = err

# print_table(1,95,1/32,10,1,rel=true)
def print_table(j0,j1,h,d,i,rel=false):
   ncols = d+1 + i+1
   print ("static const double C[" + str(j1-j0) + "][" + str(ncols) + "] = {")
   h = RR(h)
   for j in range(j0,j1):
      p, err = sollya_approx(j*2*h,h,d,i,rel=rel)
      s = "   {"
      for k in [0..i]: # double-double coefficients
         c = p[k]
         hi = RR(c)
         lo = RR(c-c.parent()(hi))
         s += get_hex(hi) + ", " + get_hex(lo) + ", "
      for k in [i+1..d]: # double coefficients
         c = p[k]
         s += get_hex(c)
         if k<d:
            s += ", "
      s += "}, /* j=" + str(j) + " " + str(float(err)) + " */"
      print (s)
   print ("};")

def extract_coeffs(p,d,i):
   l = []
   for k in [0..i]:
      c = p[k]
      hi = RR(c)
      l.append(hi)
      lo = RR(c-c.parent()(hi))
      l.append(lo)
   for k in [i+1..d]:
      l.append(p[k])
   return l

# return the 'ulp' of the interval x, i.e., max(ulp(t)) for t in x
# this internal routine is used below
def RIFulp(x):
   return max(x.lower().ulp(),x.upper().ulp())

# analyse the maximal rounding error for subrange j
# h=1/32 is hard-coded, same for d=10 and i=1
# analyze_p(1)
# 5.3078398332516484608512456856802479429227947567592236369060e-19
# analyze_p(94)
# 9.8035501070559374685251505050884232636461622258809009987654e-32
def analyze_p(j,zmin=-1/32,zmax=1/32,verbose=false):
   d = 10
   i = 1
   h = RR(1/32)
   p, err0 = sollya_approx(j*2*h,h,d,i,rel=true)
   err0 = 2^-err0
   if verbose:
      print ("err0=", log(err0)/log(2.))
   c = extract_coeffs(p,10,1) # hard-coded
   z = RIF(zmin,zmax)
   # z2 = z*z, z4 = z2*z2
   z2 = z*z
   err_z2 = RIFulp(z2)
   z4 = z2*z2
   err_z4 = RIFulp(z4)+2*z2.abs().upper()*err_z2
   # c9 = __builtin_fma (c[12], z, c[11])
   c9 = c[12]*z+c[11]
   err1 = RIFulp(c9)*z.abs().upper()^9
   if verbose:
      print ("err1=", log(err1)/log(2.))
   # c7 = __builtin_fma (c[10], z, c[9])
   c7 = c[10]*z+c[9]
   err2 = RIFulp(c7)*z.abs().upper()^7
   if verbose:
      print ("err2=", log(err2)/log(2.))
   # c5 = __builtin_fma (c[8], z, c[7])
   c5 = c[8]*z+c[7]
   err3 = RIFulp(c5)*z.abs().upper()^5
   if verbose:
      print ("err3=", log(err3)/log(2.))
   # c3 = __builtin_fma (c[6], z, c[5])
   c3 = c[6]*z+c[5]
   err4 = RIFulp(c3)*z.abs().upper()^3
   if verbose:
      print ("err4=", log(err4)/log(2.))
   # c7 = __builtin_fma (c9, z2, c7)
   c7 = c9*z2+c7
   err5 = (RIFulp(c7)+c9.abs().upper()*err_z2)*z.abs().upper()^7
   if verbose:
      print ("err5=", log(err5)/log(2.))
   # c3 = __builtin_fma (c5, z2, c3)
   c3 = c5*z2+c3
   err6 = (RIFulp(c3)+c5.abs().upper()*err_z2)*z.abs().upper()^3
   if verbose:
      print ("err6=", log(err6)/log(2.))
   # c3 = __builtin_fma (c7, z4, c3)
   c3 = c7*z4+c3
   err7 = (RIFulp(c3)+c7.abs().upper()*err_z4)*z.abs().upper()^3
   if verbose:
      print ("err7=", log(err7)/log(2.))
   # c2 = __builtin_fma (c3, z, c[4])
   c2 = c3*z+c[4]
   err8 = RIFulp(c2)*z.abs().upper()^2
   if verbose:
      print ("err8=", log(err8)/log(2.))
   # a_mul (&th, &tl, z, c2)
   th = z*c2
   u = RIFulp(th)
   tl = RIF(-u,u)
   # fast_two_sum (h, l, c[2], th)
   h = c[2]+th
   u = RIFulp(h)
   l = RIF(-u,u)
   err9 = h.abs().upper()*2^-105*z.abs().upper()
   if verbose:
      print ("err9=", log(err9)/log(2.))
   # l += tl + c[3]
   tl += c[3]
   l += tl
   err10 = (RIFulp(tl)+RIFulp(l))*z.abs().upper()
   if verbose:
      print ("err10=", log(err10)/log(2.))
   # a_mul (&th, &tl, z, h)
   th = z*h
   u = RIFulp(th)
   tl = RIF(-u,u)
   # tl = __builtin_fma (z, l, tl)
   tl += z*l
   err11 = RIFulp(tl)
   if verbose:
      print ("err11=", log(err11)/log(2.))
   # fast_two_sum (h, l, c[0], th)
   h = c[0]+th
   u = RIFulp(h)
   l = RIF(-u,u)
   err12 = h.abs().upper()*2^-105
   if verbose:
      print ("err12=", log(err12)/log(2.))
   # l += tl + c[1]
   tl += c[1]
   l += tl
   err13 = RIFulp(tl)+RIFulp(l)
   if verbose:
      print ("err13=", log(err13)/log(2.))
   err = err1+err2+err3+err4+err5+err6+err7+err8+err9+err10+err11+err12+err13
   if verbose:
      print ("err=", log(err)/log(2.))
   # convert into relative error
   if verbose:
      print ("h+l:", (h+l).abs().lower(), (h+l).abs().upper())
   err = err/(h+l).abs().lower()
   if verbose:
      print ("relerr=", log(err)/log(2.))
   # take into account err0
   R = RealField(200)
   err = (1+R(err))*(1+R(err0))-1
   return err

# call analyze_p() with 2^k subintervals
# Analyze_p(1,k=8)
# -1/32 -127/4096 3.3821100332499256013659495108510157401457307061110199316943e-19
def Analyze_p(j,k=0):
   maxerr = 0
   h = 1/16/2^k
   for i in range(2^k):
      zmin = -1/32+i*h
      zmax = zmin+h
      err = analyze_p(j,zmin=zmin,zmax=zmax)
      if err>maxerr:
         maxerr = err
         print (zmin, zmax, err)
   return maxerr

# same as analyze_p(), for special code for 0 <= z < 1/16
# analyze0()
def analyze0(zmin=0,zmax=1/16,verbose=false):
   err0 = 2^-68.935
   z = RIF(zmin,zmax)
   c0 = ["0x1.20dd750429b6dp+0", "0x1.1ae3a7862d9c4p-56", "-0x1.812746b0379e7p-2", "0x1.f1a64d72722a2p-57", "0x1.ce2f21a042b7fp-4", "-0x1.b82ce31189904p-6", "0x1.565bbf8a0fe0bp-8", "-0x1.bf9f8d2c202e4p-11"]
   c0 = [RR(x,16) for x in c0]
   # z2 = z*z
   z2 = z*z
   err_z2 = RIFulp(z2)
   # z4 = z2*z2
   z4 = z2*z2
   err_z4 = RIFulp(z4)+2*z2.abs().upper()*err_z2
   # c9 = __builtin_fma (c0[7], z2, c0[6])
   c9 = c0[7]*z2+c0[6]
   err1 = (RIFulp(c9)+c0[7]*err_z2)*z.abs().upper()^9
   if verbose:
      print ("err1=", log(err1)/log(2.))
   # c5 = __builtin_fma (c0[5], z2, c0[4])
   c5 = c0[5]*z2+c0[4]
   err2 = (RIFulp(c5)+c0[5]*err_z2)*z.abs().upper()^5
   if verbose:
      print ("err2=", log(err2)/log(2.))
   # c5 = __builtin_fma (c9, z4, c5)
   c5 = c9*z4+c5
   err3 = (RIFulp(c5)+c9.abs().upper()*err_z4)*z.abs().upper()^5
   if verbose:
      print ("err3=", log(err3)/log(2.))
   # a_mul (&th, &tl, z2, c5)
   th = z2*c5
   u = RIFulp(th)
   tl = RIF(-u,u)
   # fast_two_sum (h, l, c0[2], th)
   h = c0[2]+th
   u = RIFulp(h)
   l = RIF(-u,u)
   err4 = h.abs().upper()*2^-105*z.abs().upper()^3
   if verbose:
      print ("err4=", log(err4)/log(2.))
   # l += tl + c0[3]
   tl += c0[3]
   l += tl
   err5 = (RIFulp(tl)+RIFulp(l))*z.abs().upper()^3
   if verbose:
      print ("err5=", log(err5)/log(2.))
   # a_mul (&th, &tl, z2, h)
   th = z2*h
   u = RIFulp(th)
   tl = RIF(-u,u)
   # tl = __builtin_fma (z2, l, tl)
   tl = z2*l+tl
   err6 = (RIFulp(tl)+err_z2*l.abs().upper())*z.abs().upper()
   if verbose:
      print ("err6=", log(err6)/log(2.))
   # fast_two_sum (h, l, c0[0], th)
   h = c0[0]+th
   u = RIFulp(h)
   l = RIF(-u,u)
   err7 = h.abs().upper()*2^-105*z.abs().upper()
   if verbose:
      print ("err7=", log(err7)/log(2.))
   # l += tl + c0[1]
   tl += c0[1]
   l += tl
   err8 = (RIFulp(tl)+RIFulp(l))*z.abs().upper()
   if verbose:
      print ("err8=", log(err8)/log(2.))
   # a_mul (h, &tl, h, z)
   h = h*z
   u = RIFulp(h)
   tl = RIF(-u,u)
   # *l = __builtin_fma (*l, z, tl)
   l = l*z+tl
   err9 = RIFulp(l)
   if verbose:
      print ("err9=", log(err9)/log(2.))
   err = err1+err2+err3+err4+err5+err6+err7+err8+err9
   # convert into relative error
   err = err/(h+l).abs().lower()
   # take into account err0
   R = RealField(200)
   err = (1+R(err))*(1+R(err0))-1
   return err



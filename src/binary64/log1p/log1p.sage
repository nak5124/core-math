# check if log1p(x) rounds identically to log(x) for all rounding modes
def compare_log(x,verbose=false):
   X = x.exact_rational()
   prec = 200
   for r in ['RNDN','RNDZ','RNDU','RNDD']:
      R = RealField(53,rnd=r)
      y = R(n(log(X),prec))
      z = R(n(log(1+X),prec))
      if y!=z:
         if verbose:
            print (get_hex(x),r)
         return false
   return true

def fast_two_sum (a,b):
   hi = a + b
   e = hi - a
   lo = b - e
   return hi,lo

# check if fast_two_sum(x,1) is exact when |x| < 1
# fails with 2^-54 and rounding to nearest
def check_fast_two_sum(e):
   while true:
      x = RR.random_element()*2^e
      hi, lo = fast_two_sum(x,1)
      a = 1 + x.exact_rational()
      b = hi.exact_rational() + lo.exact_rational()
      if a!=b:
         print ("x=", get_hex(x))
         raise ValueError

# return the 'ulp' of the interval x, i.e., max(ulp(t)) for t in x
# this internal routine is used below
def RIFulp(x):
   return max(x.lower().ulp(),x.upper().ulp())

# compute the maximal absolute rounding error for p_1() over [zmin,zmax]
# if rel=true, takes into account the relative error for the Sollya polynomial
# analyze_p1(RR(2^-1000),RR(2^-999),rel=true,verbose=true)
# (1.14506750044579e-322, -1.?e-602)
# Xmax=0.00212097167968735
# analyze_p1(Xmax/2,Xmax,rel=true,verbose=true)
# (2.76450182935071e-24, 1.?e-6)
def analyze_p1(zmin,zmax,rel=true,verbose=false):
   P = ["0x1p0","-0x1p-1","0x1.555555555555p-2","-0x1.fffffffff572dp-3","0x1.999999a2d7868p-3","-0x1.5555c0d31b08ep-3","0x1.2476b9058e396p-3"]
   P = [RR(x,16) for x in P]
   z = RIF(zmin,zmax)
   if rel==false:
      err0 = 2^-79.592 # absolute error from Sollya polynomial
   else:
      err0 = 2^-70.467 # relative error from Sollya polynomial
      # for |x| < 0.00212097167968735, we have |log(1+x)/x| < 1.0011
      err0 = err0*1.0011*z.abs().upper()
   if verbose:
      print ("err0=", log(err0)/log(2.))
   # a_mul (&z2h, &z2l, z, z)
   z2h = z*z
   u = RIFulp(z2h)
   z2l = RIF(-u,u)
   err_z2h = RIFulp(z2h)
   # p56 = __builtin_fma (P[6], z, P[5])
   p56 = P[6]*z+P[5]
   err1 = RIFulp(p56)*z.abs().upper()^6
   if verbose:
      print ("err1=", log(err1)/log(2.))
   # p34 = __builtin_fma (P[4], z, P[3])
   p34 = P[4]*z+P[3]
   err2 = RIFulp(p34)*z.abs().upper()^4
   if verbose:
      print ("err2=", log(err2)/log(2.))
   # ph = __builtin_fma (p56, z2h, p34)
   ph = p56*z2h+p34
   err3 = (RIFulp(ph)+p56.abs().upper()*err_z2h)*z.abs().upper()^4
   if verbose:
      print ("err3=", log(err3)/log(2.))
   # ph = __builtin_fma (ph, z, P[2])
   ph_in = ph
   ph = ph*z+P[2]
   err4 = RIFulp(ph)*z.abs().upper()*z.abs().upper()^3
   if verbose:
      print ("err4=", log(err4)/log(2.))
   # ph *= z2h
   ph_in = ph
   ph = ph*z2h
   err5 = (RIFulp(ph)+ph_in.abs().upper()*err_z2h)*z.abs().upper()
   if verbose:
      print ("err5=", log(err5)/log(2.))
   # fast_two_sum (h, l, -0.5 * z2h, ph * z)
   t = ph*z
   h = RIF(-0.5)*z2h+t
   u = RIFulp(h)
   l = RIF(-u,u)
   err6 = RIFulp(t)+h.abs().upper()*2^-105   
   if verbose:
      print ("err6=", log(err6)/log(2.))
   # *l += -0.5 * z2l
   l = l+RIF(-0.5)*z2l
   err7 = RIFulp(l)
   if verbose:
      print ("err7=", log(err7)/log(2.))
   err = err0+err1+err2+err3+err4+err5+err6+err7
   if verbose:
      print ("err=", log(err)/log(2.))
   return err, h, l

# analyse the error of cr_log1p_fast() for |x| < 0.00212097167968735      
# analyze_x_plus_p1()
# e= 0 err= -78.3681972481953
# analyze_x_plus_p1(rel=true)
# e= 0 err= -68.3751331997836
def analyze_x_plus_p1(rel=false,verbose_e=[],Xmax=0.00212097167968735):
   maxerr = 0
   assert 0<Xmax<=1, "0<Xmax<=1"
   for e in range(0,-1075,-1):
      # consider 2^e*[Xmax/2,Xmax]
      xmin = Xmax*2^(e-1)
      xmax = Xmax*2^e
      if xmax<2^-1074:
         break
      # p_1 (h, &lo, x)
      if e in verbose_e:
         err1, h, lo = analyze_p1(xmin,xmax,rel=rel,verbose=true)
         print ("err1=", log(err1)/log(2.))
      else:
         err1, h, lo = analyze_p1(xmin,xmax,rel=rel)
      x = RIF(xmin,xmax)
      # fast_two_sum (h, l, x, *h)
      h = x + h
      u = RIFulp(h)
      l = RIF(-u,u)
      err2 = h.abs().upper()*2^-105
      if e in verbose_e:
         print ("err2=", log(err2)/log(2.))
      # *l += lo
      l += lo
      err3 = RIFulp(l)
      err = err1+err2+err3
      if rel:
         err = err/h.abs().lower()
      if err>maxerr:
         print ("e=", e, "err=", log(err)/log(2.))
         maxerr = err
      # same for negative values
      xmin,xmax = -xmax,-xmin
      # p_1 (h, &lo, x)
      if e in verbose_e:
         err1, h, lo = analyze_p1(xmin,xmax,rel=rel,verbose=true)
         print ("err1=", log(err1)/log(2.))
      else:
         err1, h, lo = analyze_p1(xmin,xmax,rel=rel)
      x = RIF(xmin,xmax)
      # fast_two_sum (h, l, x, *h)
      h = x + h
      u = RIFulp(h)
      l = RIF(-u,u)
      err2 = h.abs().upper()*2^-105
      if e in verbose_e:
         print ("err2=", log(err2)/log(2.))
      # *l += lo
      l += lo
      err3 = RIFulp(l)
      err = err1+err2+err3
      if rel:
         err = err/h.abs().lower()
      if err>maxerr:
         print ("e=", e, "err=", log(err)/log(2.))
         maxerr = err

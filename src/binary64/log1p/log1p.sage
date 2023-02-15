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
# (1.34747385904737e-21, -0.00001?)
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
   # z2 = z * z
   z2 = z*z
   err_z2 = RIFulp(z2)
   # p56 = __builtin_fma (P[6], z, P[5])
   p56 = P[6]*z+P[5]
   err1 = RIFulp(p56)*z2.abs().upper()^3
   if verbose:
      print ("err1=", log(err1)/log(2.))
   # p34 = __builtin_fma (P[4], z, P[3])
   p34 = P[4]*z+P[3]
   err2 = RIFulp(p34)*z2.abs().upper()^2
   if verbose:
      print ("err2=", log(err2)/log(2.))
   # p12 = __builtin_fma (P[2], z, P[1])
   p12 = P[2]*z+P[1]
   err3 = RIFulp(p12)*z2.abs().upper()
   if verbose:
      print ("err3=", log(err3)/log(2.))
   # ph = __builtin_fma (p56, z2, p34)
   ph = p56*z2+p34
   err4 = (RIFulp(ph)+p56.abs().upper()*err_z2)*z2.abs().upper()^2
   if verbose:
      print ("err4=", log(err4)/log(2.))
   # ph = __builtin_fma (ph, z2, p12)
   ph_in = ph
   ph = ph*z2+p12
   err5 = (RIFulp(ph)+ph_in.abs().upper()*err_z2)*z2.abs().upper()
   if verbose:
      print ("err5=", log(err5)/log(2.))
   # ph *= z2
   ph_in = ph
   ph = ph*z2
   err6 = RIFulp(ph)+ph_in.abs().upper()*err_z2
   if verbose:
      print ("err6=", log(err6)/log(2.))
   err = err0+err1+err2+err3+err4+err5+err6
   if verbose:
      print ("err=", log(err)/log(2.))
   return err, ph

# analyse the error of cr_log1p_fast() for |x| < 0.00212097167968735      
# analyze_x_plus_p1()
# e= 0 err= -68.8748330313666
# analyze_x_plus_p1(rel=true)
# e= 0 err= -58.9943804324056
def analyze_x_plus_p1(rel=false,verbose_e=[]):
   Xmax = 0.00212097167968735
   maxerr = 0
   for e in range(0,-1067,-1):
      # consider 2^e*[Xmax/2,Xmax]
      xmin = Xmax*2^(e-1)
      xmax = Xmax*2^e
      if e in verbose_e:
         err1, p = analyze_p1(xmin,xmax,rel=rel,verbose=true)
         print ("err1=", log(err1)/log(2.))
      else:
         err1, p = analyze_p1(xmin,xmax,rel=rel)
      x = RIF(xmin,xmax)
      # fast_two_sum (h, l, x, p)
      h = x + p
      err2 = h.abs().upper()*2^-105
      if e in verbose_e:
         print ("err2=", log(err2)/log(2.))
      err = err1+err2
      if rel:
         err = err/h.abs().lower()
      if err>maxerr:
         print ("e=", e, "err=", log(err)/log(2.))
         maxerr = err
      xmin,xmax = -xmax,-xmin
      if e in verbose_e:
         err1, p = analyze_p1(xmin,xmax,rel=rel,verbose=true)
         print ("err1=", log(err1)/log(2.))
      else:
         err1, p = analyze_p1(xmin,xmax,rel=rel)
      x = RIF(xmin,xmax)
      # fast_two_sum (h, l, x, p)
      h = x + p
      err2 = h.abs().upper()*2^-105
      if e in verbose_e:
         print ("err2=", log(err2)/log(2.))
      err = err1+err2
      if rel:
         err = err/h.abs().lower()
      if err>maxerr:
         print ("e=", e, "err=", log(err)/log(2.))
         maxerr = err

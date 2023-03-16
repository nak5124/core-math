def a_mul(a,b):
   hi = a*b
   lo = fma(a,b,-hi)
   return hi, lo

def d_mul(ah,al,bh,bl):
   hi, s = a_mul (ah, bh)
   t = fma (al, bh, s)
   lo = fma (ah, bl, t)
   return hi, lo

def d_inv(bh,bl):
   R = bh.parent()
   hi = R(1) / bh
   e = fma (-hi, bh, R(1))
   e = fma (-hi, bl, e)
   lo = hi*e
   return hi, lo

def d_div(ah,al,bh,bl):
   hi, lo = d_inv (bh, bl)
   return d_mul (ah, al, hi, lo)

# ah= 0x1.014a3e61bc114p-1 al= 0x1.f56cd49187432p-54 bh= 0x1.0ddd63e28aadep-1 bl= -0x1.f7fdbb36e232cp-54 err= -102.052325337202
def check_d_div(K=10^6,rnd='RNDN'):
   maxerr = 0
   R = RealField(53,rnd=rnd)
   for k in range(K):
      ah = R.random_element()
      al = R.random_element()*ah.ulp()
      bh = R.random_element()
      bl = R.random_element()*bh.ulp()
      hi, lo = d_div(ah,al,bh,bl)
      a = ah.exact_rational()+al.exact_rational()
      b = bh.exact_rational()+bl.exact_rational()
      r = a/b
      err = n((hi.exact_rational()+lo.exact_rational())/r-1,200)
      err = abs(err)
      if err>maxerr:
         print ("ah=", get_hex(ah), "al=", get_hex(al), "bh=", get_hex(bh), "bl=", get_hex(bl), "err=", log(err)/log(2.))
         maxerr = err
         
# return the 'ulp' of the interval x, i.e., max(ulp(t)) for t in x
# this internal routine is used below
def RIFulp(x):
   return max(x.lower().ulp(),x.upper().ulp())

def analyze_eval_T(wmin,wmax,rel=false,verbose=false):
   w = RIF(wmin,wmax)
   if rel==false:
      err0 = 2^-80.528
   else: # |tanh(x)/x| < 1
      err0 = 2^-71.98*w.abs().upper()
   T = ["0x1p0","-0x1.5555555555553p-2","0x1.111111103f43cp-3","-0x1.ba18e77264096p-5"]
   T = [RR(x,16) for x in T]
   # z = w * w
   z = w*w
   errz = RIFulp(z)
   # *h = __builtin_fma (TT[3], z, TT[2])
   h = T[3]*z+T[2]
   err1 = (RIFulp(h)+T[3]*errz)*w.abs().upper()^5
   if verbose:
      print ("err1=", log(err1)/log(2.))
   # *h = __builtin_fma (*h, z, TT[1])
   hin = h
   h = hin*z+T[1]
   err2 = (RIFulp(h)+hin.abs().upper()*errz)*w.abs().upper()^3
   if verbose:
      print ("err2=", log(err2)/log(2.))
   # *h = *h * z
   hin = h
   h = h*z
   err3 = (RIFulp(h)+hin.abs().upper()*errz)*w.abs().upper()
   if verbose:
      print ("err3=", log(err3)/log(2.))
   # fast_two_sum (h, l, w, *h * w)
   h = w + h*w
   err4 = 2^-105*h.abs().upper()
   if verbose:
      print ("err4=", log(err4)/log(2.))
   err = err0+err1+err2+err3+err4
   if rel:
      err = err/h.abs().lower()
   if verbose:
      print ("err=", log(err)/log(2.))
   return err

# analyze_eval_T_all()
# 2.47926464106338e-23
# analyze_eval_T_all(rel=true)
# 9.35004837242004e-21
def analyze_eval_T_all(rel=false):
   maxerr = 0
   e = 0
   while true:
      wmax = 0.00543*2^e
      wmin = wmax/2
      err = analyze_eval_T(wmin,wmax,rel=rel)
      if err>maxerr:
         print ("e=", e, "err=", err)
         maxerr = err
      err = analyze_eval_T(-wmax,-wmin,rel=rel)
      if err>maxerr:
         print ("e=", e, "err=", err)
         maxerr = err
      if wmin<2^-1074:
         break
      e = e-1
   return maxerr

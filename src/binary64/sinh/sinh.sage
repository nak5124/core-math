def sinh_ok(x,k):
   R = RealField(53+k)
   y = sinh(R(x))
   return y.exact_rational()==RR(y).exact_rational()

def cosh_ok(x,k):
   R = RealField(53+k)
   y = cosh(R(x))
   return y.exact_rational()==RR(y).exact_rational()

def exp_ok(x,k):
   R = RealField(53+k)
   y = cosh(R(x))
   return y.exact_rational()==RR(y).exact_rational()

# compute (xi,sinh(xi),exp(-xi)) for xi near i*2^8/0x1.70f77fc88ae3cp6
# such that sinh(xi) and cosh(xi)=sinh(xi)+exp(-xi) are accurate to 53+k bits
# build_table_T(k=16)
# 2.35540937865153e-8
def build_table_T(k=16):
   magic = RR("0x1.70f77fc88ae3cp6",16)
   h = 1/magic
   maxerr = 0
   print ("static const double T[256][3] = {")
   for i in range(256):
      if i==0:
         x0 = RR(0)
      else:
         x0 = x1 = i*2^8*h
         # compute the exponent difference between sinh(x0) and exp(-x0)
         r = sinh(x0)/exp(-x0)
         diff = ZZ(ceil(log(r)/log(2.)))
         # we need exp(-xi) accurate to only k-diff bits
         k1 = max(0,k-diff)
         while true:
            if sinh_ok(x0,k) and exp_ok(-x0,k1):
               break
            x0 = x0.nextbelow()
            x1 = x1.nextabove()
            if sinh_ok(x1,k) and exp_ok(-x1,k1):
               x0 = x1
               break
      maxerr = max(maxerr,abs(x0-i*2^8*h))
      print ("   {" + get_hex(x0) + ", " + get_hex(sinh(x0)) + ", " + get_hex(exp(-x0)) + "}, /* i=" + str(i) + " */")
   print ("};")
   return maxerr

# compute (xi,sinh(xi),cosh(xi)) for xi near i/0x1.70f77fc88ae3cp6
# such that sinh(xi) and cosh(xi) are accurate to 53+k bits
# return the maximal difference |xi - i/0x1.70f77fc88ae3cp6|
# build_table_U(k=11)
# 3.74086206456070e-9
# build_table_U(k=12)
# 1.91293030304962e-8
def build_table_U(k=10):
   magic = RR("0x1.70f77fc88ae3cp6",16)
   h = 1/magic
   maxerr = 0
   print ("static const double U[256][3] = {")
   for i in range(0,256):
      if i==0:
         x0 = RR(0)
      else:
         x0 = x1 = i*h
         while true:
            if sinh_ok(x0,k) and cosh_ok(x0,k):
               break
            x0 = x0.nextbelow()
            x1 = x1.nextabove()
            if sinh_ok(x1,k) and cosh_ok(x1,k):
               x0 = x1
               break
      maxerr = max(maxerr,abs(x0-i*h))
      print ("   {" + get_hex(x0) + ", " + get_hex(sinh(x0)) + ", " + get_hex(cosh(x0)) + "}, /* i=" + str(i) + " */")
   print ("};")
   return maxerr

def build_table_Tl(T0,T1,T2):
   print ("static const double Tl[256][2] = {")
   for i in range(0,256):
      x = RR(T0[i],16)
      s = RR(T1[i],16)
      e = RR(T2[i],16)
      X = x.exact_rational()
      ls = RR(n(sinh(X)-s.exact_rational(),200))
      le = RR(n(cosh(X)-sinh(X)-e.exact_rational(),200))
      print ("   {" + get_hex(ls) + ", " + get_hex(le) + "}, /* i=" + str(i) + " */")
   print ("};")
            
def build_table_Ul(U0,U1,U2):
   print ("static const double Ul[256][2] = {")
   for i in range(0,256):
      x = RR(U0[i],16)
      s = RR(U1[i],16)
      c = RR(U2[i],16)
      X = x.exact_rational()
      ls = RR(n(sinh(X)-s.exact_rational(),200))
      lc = RR(n(cosh(X)-c.exact_rational(),200))
      print ("   {" + get_hex(ls) + ", " + get_hex(lc) + "}, /* i=" + str(i) + " */")
   print ("};")
            
# return the 'ulp' of the interval x, i.e., max(ulp(t)) for t in x
# this internal routine is used below
def RIFulp(x):
   return max(x.lower().ulp(),x.upper().ulp())

def analyze_eval_S(wmin,wmax,rel=false,verbose=false):
   w = RIF(wmin,wmax)
   if rel==false:
      err0 = 2^-83.263
   else: # |sinh(x)/x| < 1.000005 for |x| < 0.00543
      err0 = 2^-74.818*1.000005*w.abs().upper()
   S = ["0x1p0","0x1.5555555555555p-3","0x1.11111111869d4p-7","0x1.a01061b363a81p-13"]
   S = [RR(x,16) for x in S]
   # z = w * w
   z = w*w
   errz = RIFulp(z)
   # *h = __builtin_fma (S[3], z, S[2])
   h = S[3]*z+S[2]
   err1 = (RIFulp(h)+S[3]*errz)*w.abs().upper()^5
   if verbose:
      print ("err1=", log(err1)/log(2.))
   # *h = __builtin_fma (*h, z, S[1])
   hin = h
   h = hin*z+S[1]
   err2 = (RIFulp(h)+hin.abs().upper()*errz)*w.abs().upper()^3
   if verbose:
      print ("err2=", log(err2)/log(2.))
   # *h = *h * z
   h = h*z
   err3 = RIFulp(h)*w.abs().upper()
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

# analyze_eval_S_all()
# 9.12931943224523e-24
# analyze_eval_S_all(rel=true)
# 3.39086502646348e-21
def analyze_eval_S_all(rel=false):
   maxerr = 0
   e = 0
   while true:
      wmax = 0.00543*2^e
      wmin = wmax/2
      err = analyze_eval_S(wmin,wmax,rel=rel)
      if err>maxerr:
         print ("e=", e, "err=", err)
         maxerr = err
      err = analyze_eval_S(-wmax,-wmin,rel=rel)
      if err>maxerr:
         print ("e=", e, "err=", err)
         maxerr = err
      if wmin<2^-1074:
         break
      e = e-1
   return maxerr

# analyze_eval_C(verbose=true)
# err1= -87.0993273681392
# err2= -68.0496623387826
# err3= -104.415036317674
# err= -68.0496596801630
# 3.27386391283508e-21
def analyze_eval_C(verbose=false):
   err0 = 2^-81.152
   w = RIF(-0.00543,0.00543)
   C = ["0x1p0","0x1p-1","0x1.5555555554e2ep-5","0x1.6c16d52a52a35p-10"]
   C = [RR(x,16) for x in C]
   # z = w * w
   z = w*w
   errz = RIFulp(z)
   # *h = __builtin_fma (C[3], z, C[2])
   h = C[3]*z+C[2]
   err1 = (RIFulp(h)+C[3]*errz)*w.abs().upper()^4
   if verbose:
      print ("err1=", log(err1)/log(2.))
   # *h = __builtin_fma (*h, z, C[1])
   hin = h
   h = hin*z+C[1]
   err2 = (RIFulp(h)+hin.abs().upper()*errz)*w.abs().upper()^2
   if verbose:
      print ("err2=", log(err2)/log(2.))
   # fast_two_sum (h, l, 1.0, *h)
   h = RR(1) + h
   err3 = 2^-105*h.abs().upper()
   if verbose:
      print ("err3=", log(err3)/log(2.))
   err = err0+err1+err2+err3
   if verbose:
      print ("err=", log(err)/log(2.))
   return err
   
def check_sw(w,swh,swl):
   W = w.exact_rational()
   H = swh.exact_rational()
   L = swl.exact_rational()
   e = n(sinh(W)-H-L,200)
   return log(abs(e/swh))/log(2.)

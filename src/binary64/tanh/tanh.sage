def fast_two_sum(a,b):
   hi = a+b
   e = hi-a
   lo = b-e
   return hi, lo

def two_sum(a,b):
   s = a+b
   ap = s-b
   bp = s-ap
   da = a-ap
   db = b-bp
   t = da+db
   return s,t

def a_mul(a,b):
   hi = a*b
   lo = fma(a,b,-hi)
   return hi, lo

def d_mul(ah,al,bh,bl):
   hi, s = a_mul (ah, bh)
   t = fma (al, bh, s)
   lo = fma (ah, bl, t)
   return hi, lo

def d_mul_acc2(ah,al,bh,bl):
   hi, cl1 = a_mul(ah,bh)
   tl = ah * bl
   cl2 = fma(al, bh, tl)
   cl3 = cl1 + cl2
   return fast_two_sum (hi, cl3)

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

def s_mul_acc1(a,bh,bl):
   hi, cl1 = a_mul (a, bh)
   cl2 = a * bl
   hi, lo = fast_two_sum (hi, cl2)
   tl2 = lo + cl1
   return fast_two_sum (hi, tl2)

def s_mul_acc3(a,bh,bl):
   hi, cl1 = a_mul(a,bh)
   cl3 = fma(a,bl,cl1)
   return fast_two_sum (hi, cl3)

def d_div_acc1(xh,xl,yh,yl):
   th = xh / yh
   rh, rl = s_mul_acc1 (th, yh, yl)
   pih, pil = two_sum (xh, -rh)
   delta_h = pil - rl
   delta_l = delta_h + xl
   delta = pih + delta_l
   tl = delta / yh
   return fast_two_sum (th, tl)

def d_div_acc2(xh,xl,yh,yl):
   th = xh / yh
   rh, rl = s_mul_acc1 (th, yh, yl)
   pih = xh - rh
   delta_l = xl - rl
   delta = pih + delta_l
   tl = delta / yh
   return fast_two_sum (th, tl)

def fast_sum_acc(a,bh,bl):
   sh, sl = two_sum(bh,a)
   v = bl+sl
   return fast_two_sum (sh, v)

def d_div_acc3(xh,xl,yh,yl):
   R = xh.parent()
   th = R(1)/yh
   rh = fma(-yh,th,R(1))
   rl = -yl * th
   eh, el = fast_two_sum (rh, rl)
   delta_h, delta_l = s_mul_acc3 (th, eh, el)
   mh, ml = fast_sum_acc (th, delta_h, delta_l)
   return d_mul_acc2 (xh, xl, mh, ml)

# check_d_div(K=10^6,rnd='RNDN')
# ah= 0x1.014a3e61bc114p-1 al= 0x1.f56cd49187432p-54 bh= 0x1.0ddd63e28aadep-1 bl= -0x1.f7fdbb36e232cp-54 err= 15.4300909741780
# avg. err= 1.14459486977720
# check_d_div(K=10^6,rnd='RNDZ')
# ah= -0x1.00933ddb51e1cp-1 al= -0x1.eb3917eafca2ap-54 bh= -0x1.426a3514adee4p-1 bl= 0x1.fe34e12a06732p-54 err= 30.6708140509637
# avg. err= 4.33318176146531
# check_d_div(K=10^6,rnd='RNDN',algo=d_div_acc1)
# ah= 0x1.2f7e7879b5f74p-1 al= 0x1.af8bc1f1bcf0ep-54 bh= 0x1.0bff358d60384p-1 bl= -0x1.f5a2c6a9b41bep-54 err= 15.2999565566288
# avg. err= 0.986063329074064
# check_d_div(K=10^6,rnd='RNDN',algo=d_div_acc2)
# ah= 0x1.22aa4690bdd0ap-1 al= 0x1.6878cca074226p-54 bh= 0x1.09b548efe132p-2 bl= -0x1.f33ded8ee1422p-55 err= 14.2802591373840
# avg. err= 0.987773577175927
# check_d_div(K=10^6,rnd='RNDN',algo=d_div_acc3)
# ah= 0x1.3cf7a1a1cfbb8p-1 al= -0x1.78eeeffc6ef32p-54 bh= -0x1.05baa81804abp-2 bl= -0x1.f9b471a10ff76p-55 err= 8.95713681827043
# avg. err= 1.08553630618232
def check_d_div(K=10^6,rnd='RNDN',algo=d_div):
   maxerr = 0
   avgerr = 0
   R = RealField(53,rnd=rnd)
   u = 2^-53.
   for k in range(K):
      ah = R.random_element()
      al = R.random_element()*ah.ulp()
      bh = R.random_element()
      bl = R.random_element()*bh.ulp()
      hi, lo = algo(ah,al,bh,bl)
      a = ah.exact_rational()+al.exact_rational()
      b = bh.exact_rational()+bl.exact_rational()
      r = a/b
      err = n((hi.exact_rational()+lo.exact_rational())/r-1,200)
      err = abs(err)/u^2
      avgerr += err
      if err>maxerr:
         print ("ah=", get_hex(ah), "al=", get_hex(al), "bh=", get_hex(bh), "bl=", get_hex(bl), "err=", err)
         maxerr = err
   avgerr /= K
   print ("avg. err=", avgerr)
         
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
# [((8183583624766079, -79), (5365348628892106, -48))]
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

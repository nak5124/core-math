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
# ah= -0x1.0a7f842d1d77ep-1 al= 0x1.56ac275541066p-55 bh= 0x1.1450e2df52752p-1 bl= 0x1.d2c2cc7488ea8p-55 err= 6.69024198067564
# avg. err= 0.560541754000017
# check_d_div(K=10^6,rnd='RNDN',worst=true)
# ah= 0x1.13ef3bfc034ecp-1 al= 0x1p-54 bh= -0x1.a1e103c44907cp-1 bl= 0x1p-54 err= 8.71261829290411
# check_d_div(K=10^6,rnd='RNDZ')
# ah= 0x1.11893f386fap-7 al= 0x1.aa0e8e2e3c4bcp-61 bh= -0x1.e73756bc4ef46p-1 bl= 0x1.f5cf9a97b4b7cp-55 err= 24.7338189683173
# avg. err= 3.64249951725914
# check_d_div(K=10^6,rnd='RNDN',algo=d_div_acc1) # claimed 15u^2
# ah= -0x1.0ecc6af81f4acp-1 al= -0x1.91e410c63551ap-55 bh= 0x1.022185e49a19p-3 bl= -0x1.beb9dc847335ap-57 err= 6.51290290645275
# avg. err= 0.418192527364857
# check_d_div(K=10^6,rnd='RNDN',algo=d_div_acc1,worst=true)
# ah= -0x1.0d635e390f03p-1 al= 0x1p-54 bh= -0x1.068e39e3db3cp-6 bl= -0x1p-59 err= 6.41238254052790
# check_d_div(K=10^6,rnd='RNDN',algo=d_div_acc2) # claimed 15u^2
# ah= -0x1.1e99e0476fd32p-1 al= 0x1.ff23473ea10e6p-55 bh= -0x1.0c6766d79d598p-3 bl= -0x1.d547f7bdd0aep-57 err= 6.26011749470428
# avg. err= 0.417522466227338
# check_d_div(K=10^6,rnd='RNDN',algo=d_div_acc2,worst=true)
# ah= 0x1.04248b8f98064p-1 al= -0x1p-54 bh= 0x1.28313836c7dcp-4 bl= 0x1p-57 err= 6.28648758693814
# check_d_div(K=10^7,rnd='RNDN',algo=d_div_acc3) # claimed 9.8u^2
# ah= 0x1.0907f3c15df88p-2 al= 0x1.b566aa7e9c0f4p-56 bh= -0x1.fa922e2e59498p-3 bl= -0x1.f89c70a395d6ep-57 err= 5.61606432576871
# avg. err= 0.549085609379064
# check_d_div(K=10^6,rnd='RNDN',algo=d_div_acc3,worst=true)
# ah= -0x1.0078454fea19ep-1 al= -0x1p-54 bh= 0x1.ed4bd283fap-9 bl= 0x1p-62 err= 5.27837842708530
def check_d_div(K=10^6,rnd='RNDN',algo=d_div,worst=false):
   maxerr = 0
   avgerr = 0
   R = RealField(53,rnd=rnd)
   u = 2^-53.
   for k in range(K):
      ah = R.random_element()
      if worst==false:
         al = R.random_element()*ah.ulp()/2
      else:
         if abs(ah.exact_rational().numer()) < 2^52:
            al = ah.ulp()/2
         else:
            al = (ah.ulp()/2).nextbelow()
         if randint(0,1)==1:
            al = -al
      bh = R.random_element()
      if worst==false:
         bl = R.random_element()*bh.ulp()/2
      else:
         if abs(bh.exact_rational().numer()) < 2^52:
            bl = bh.ulp()/2
         else:
            bl = (bh.ulp()/2).nextbelow()
         if randint(0,1)==1:
            bl = -bl
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

def print_T2():
   print ("static const long double T2[32][2] = {")
   R = RealField(64)
   maxerr = 0
   for i in range(2^5):
      x = n(2^(i/2^5), 200)
      h = R(x)
      l = R(x-h.exact_rational())
      print ("   {" + get_hex(h) + "L, " + get_hex(l) + "L},")
      maxerr = max(maxerr,abs(x-h.exact_rational()-l.exact_rational())/x)
   print ("};")
   print (log(maxerr)/log(2.))

def print_T2fast():
   print ("static const double T2fast[32][2] = {")
   R = RealField(53)
   maxerr = 0
   for i in range(2^5):
      x = n(2^(i/2^5), 200)
      h = R(x)
      l = R(x-h.exact_rational())
      print ("   {" + get_hex(h) + "L, " + get_hex(l) + "L},")
      maxerr = max(maxerr,abs(x-h.exact_rational()-l.exact_rational()))
   print ("};")
   print (log(maxerr)/log(2.))

def print_T1():
   print ("static const long double T1[32][2] = {")
   R = RealField(64)
   maxerr = 0
   for i in range(2^5):
      x = n(2^(i/2^10), 200)
      h = R(x)
      l = R(x-h.exact_rational())
      print ("   {" + get_hex(h) + "L, " + get_hex(l) + "L},")
      maxerr = max(maxerr,abs(x-h.exact_rational()-l.exact_rational())/x)
   print ("};")
   print (log(maxerr)/log(2.))

def print_T1fast():
   print ("static const double T1fast[32][2] = {")
   R = RealField(53)
   maxerr = 0
   for i in range(2^5):
      x = n(2^(i/2^10), 200)
      h = R(x)
      l = R(x-h.exact_rational())
      print ("   {" + get_hex(h) + "L, " + get_hex(l) + "L},")
      maxerr = max(maxerr,abs(x-h.exact_rational()-l.exact_rational()))
   print ("};")
   print (log(maxerr)/log(2.))

def print_T0():
   print ("static const long double T0[32][2] = {")
   R = RealField(64)
   maxerr = 0
   for i in range(2^5):
      x = n(2^(i/2^15), 200)
      h = R(x)
      l = R(x-h.exact_rational())
      print ("   {" + get_hex(h) + "L, " + get_hex(l) + "L},")
      maxerr = max(maxerr,abs(x-h.exact_rational()-l.exact_rational())/x)
   print ("};")
   print (log(maxerr)/log(2.))

def print_T0fast():
   print ("static const double T0fast[32][2] = {")
   R = RealField(53)
   maxerr = 0
   for i in range(2^5):
      x = n(2^(i/2^15), 200)
      h = R(x)
      l = R(x-h.exact_rational())
      print ("   {" + get_hex(h) + "L, " + get_hex(l) + "L},")
      maxerr = max(maxerr,abs(x-h.exact_rational()-l.exact_rational()))
   print ("};")
   print (log(maxerr)/log(2.))

# return the 'ulp' of the interval x, i.e., max(ulp(t)) for t in x
# this internal routine is used below
def RIulp(x):
   return max(x.lower().ulp(),x.upper().ulp())

def a_mul(a,b):
   hi = a*b
   u = RIulp(hi)
   lo = RIF(-u,u)
   return hi,lo

def d_mul(ah,al,bh,bl):
   hi, s = a_mul(ah,bh)
   t = al*bh+s
   lo = ah*bl+t
   return hi,lo,(RIF(RIulp(t)+RIulp(lo))+RIF(al)*RIF(bl)).abs().upper()

# given maximum absolute values, return a bound on the *absolute* error of d_mul,
# and the maximum value of lo
# err_d_mul(2.,2^-53.,2.,2^-53.)
# 4.06756404254584e-31
def err_d_mul(ahmax,almax,bhmax,blmax):
   hi = ahmax*bhmax
   s = hi.ulp()
   t = almax*bhmax+s
   lo = ahmax*blmax+t
   return t.ulp()+lo.ulp()+almax*blmax, lo

# analyze_P()
# err1= -104.999994496566
# err2= -86.9999944965657
# err3= -121.528758743543
# err4= -119.999999999999
# err5= -118.796973547466
# err6= -104.999984741210
# err7= -104.000000000001
# rel. err= -86.8876714459048
# max l= 2.22049521164608e-16
def analyze_P():
   err0 = 2^-90.627 # absolute error
   R = RealField(53)
   RI = RealIntervalField(53)
   xh = RI(-2^-16,2^-16)
   xl = RI(-2^-69,2^-69)
   p4 = RI(R("0x1.3b2a52e855b32p-7",16))
   p3 = RI(R("0x1.c6b08d7057b35p-5",16))
   p2 = RI(R("0x1.ebfbdff82c58fp-3",16))
   p1h = RI(R("0x1.62e42fefa39efp-1",16))
   p1l = RI(R("0x1.abc9c864cbd56p-56",16))
   p0 = RI(1)
   # y = p[5] * xh + p[4]
   y = p4*xh+p3
   err1 = (RIulp(p4*xh)+RIulp(y))*xh.abs().upper()^3
   print ("err1=", log(err1)/log(2.))
   # y = y * xh + p[3]
   yin = y
   y = y*xh+p2
   err2 = (RIulp(yin*xh)+RIulp(y))*xh.abs().upper()^2
   print ("err2=", log(err2)/log(2.))
   # a_mul_double (h, l, y, xh)
   h = y*xh
   u = RIulp(h)
   l = RI(-u,u)
   # fast_two_sum_double (h, &t, p[1], h)
   h = p1h+h
   u = RIulp(h)
   t = RI(-u,u)
   err3 = (h.abs().upper()*2^-105)*xh.abs().upper()
   print ("err3=", log(err3)/log(2.))
   # *l += t + p[2]
   t += p1l
   l += t
   err4 = (RIulp(t)+RIulp(l))*xh.abs().upper()
   print ("err4=", log(err4)/log(2.))
   # d_mul_double (h, l, *h, *l, xh, xl)
   h, l, err5 = d_mul(h,l,xh,xl)
   print ("err5=", log(err5)/log(2.))
   # fast_two_sum_double (h, &t, p[0], h)
   h = p0+h
   u = RIulp(h)
   t = RI(-u,u)
   err6 = h.abs().upper()*2^-105
   print ("err6=", log(err6)/log(2.))
   # *l += t
   l += t
   err7 = RIulp(l)
   print ("err7=", log(err7)/log(2.))
   # convert err0 into relative error
   err0 = err0*(h+l).abs().upper()
   err = err0+err1+err2+err3+err4+err5+err6+err7
   # convert into relative error
   err = err/(h+l).abs().lower()
   print ("rel. err=", log(err)/log(2.))
   print ("max l=", l.abs().upper())

# analyze_Pacc()
# err1= -152.999997248280
# err2= -133.999997248280
# err3= -134.000000000000
# err4= -177.415037499279
# err5= -159.433959645837
# err6= -159.415037499279
# err7= -141.570704961070
# err8= -141.415037499279
# err9= -125.415032412998
# rel. err= -125.403714690496
# max l= 1.08422046458161582e-19
def analyze_Pacc():
   err0 = 2^-133.987 # absolute error
   R = RealField(64)
   RI = RealIntervalField(64)
   x = RI(-2^-16,2^-16)
   p6 = RI(R("0x1.4309131bde9fabeap-13",16))
   p5 = RI(R("0x1.5d87fe78ad725bcep-10",16))
   p4 = RI(R("0x1.3b2ab6fba4e7729cp-7",16))
   p3h = RI(R("0x1.c6b08d704a0bf8b4p-5",16))
   p3l = RI(R("-0x1.8b4ba2fbcf44117p-70",16))
   p2h = RI(R("0x1.ebfbdff82c58ea86p-3",16))
   p2l = RI(R("0x1.e2d60dd936b9ba5ep-68",16))
   p1h = RI(R("0x1.62e42fefa39ef358p-1",16))
   p1l = RI(R("-0x1.b0e2633fe0676a9cp-67",16))
   p0 = RI(1)
   # y = p[6] * x + p[5]
   y = p6*x+p5
   err1 = (RIulp(p6*x)+RIulp(y))*x.abs().upper()^5
   print ("err1=", log(err1)/log(2.))
   # y = y * x + p[4]
   yin = y
   y = y*x+p4
   err2 = (RIulp(yin*x)+RIulp(y))*x.abs().upper()^4
   print ("err2=", log(err2)/log(2.))
   # y = y * x
   yin = y
   y = y*x
   err3a = RIulp(yin*x)
   # fast_two_sum (h, l, p3h, y)
   h = p3h+y
   u = RIulp(h)
   l = RI(-u,u)
   err3b = 2^-127*h.abs().upper()
   # l += p3l
   l += p3l
   err3c = RIulp(l)
   err3 = (err3a+err3b+err3c)*x.abs().upper()^3
   print ("err3=", log(err3)/log(2.))
   # now multiply h+l by x
   # a_mul (h, t, h, x) # exact
   h = h*x
   u = RIulp(h)
   t = RI(-u,u)
   # l = l*x+t
   lin = l
   l = l*x+t
   err4 = (RIulp(lin*x)+RIulp(l))*x.abs().upper()^2
   print ("err4=", log(err4)/log(2.))
   # add p2h+p2l
   # fast_two_sum (h, t, p2h, h)
   h = p2h+h
   u = RIulp(h)
   t = RI(-u,u)
   err5a = 2^-127*h.abs().upper()
   # l += t + p2l
   l += t+p2l
   err5 = (err5a+RIulp(t+p2l)+RIulp(l))*x.abs().upper()^2
   print ("err5=", log(err5)/log(2.))
   # now multiply h+l by x
   # a_mul (h, t, h, x) # exact
   h = h*x
   u = RIulp(h)
   t = RI(-u,u)
   # l = l*x+t
   lin = l
   l = l*x+t
   err6 = (RIulp(lin*x)+RIulp(l))*x.abs().upper()
   print ("err6=", log(err6)/log(2.))
   # add p1h+p1l
   # fast_two_sum (h, t, p1h, h)
   h = p1h+h
   u = RIulp(h)
   t = RI(-u,u)
   err7a = 2^-127*h.abs().upper()
   # l += t + p1l
   l += t+p1l
   err7 = (err7a+RIulp(t+p1l)+RIulp(l))*x.abs().upper()
   print ("err7=", log(err7)/log(2.))
   # now multiply h+l by x
   # a_mul (h, t, h, x) # exact
   h = h*x
   u = RIulp(h)
   t = RI(-u,u)
   # l = l*x+t
   lin = l
   l = l*x+t
   err8 = RIulp(lin*x)+RIulp(l)
   print ("err8=", log(err8)/log(2.))
   # add p0
   # fast_two_sum (h, t, p0, h)
   h = p0+h
   u = RIulp(h)
   t = RI(-u,u)
   err9a = 2^-127*h.abs().upper()
   # l += t
   l += t
   err9 = err9a+RIulp(l)
   print ("err9=", log(err9)/log(2.))
   # convert err0 into relative error
   err0 = err0*(h+l).abs().upper()
   err = err0+err1+err2+err3+err4+err5+err6+err7+err8+err9
   # convert into relative error
   err = err/(h+l).abs().lower()
   print ("rel. err=", log(err)/log(2.))
   print ("max l=", l.abs().upper())
   
# split binade [2^(e-1),2^e) into k blocks
def doit_bacsel(e,k,t0=None,t1=None,neg=false):
   if t0==None:
      if neg:
         t0 = -2^64+1
      else:
         t0 = 2^63
   if t1==None:
      if neg:
         t1 = -2^63+1
      else:
         t1 = 2^64
   h = ceil((t1-t0)/k)
   for i in range(k):
      u0 = t0+h*i
      u1 = min(t0+h*(i+1),t1)
      print ("./doit.sh " + str(u0) + " " + str(u1) + " 64 " + str(e) + " 64 30.5 64")

def dekker(u,v):
   R = u.parent()
   c = R("0x1.00000001p+32",16)
   up = u*c
   vp = v*c
   u1 = (u - up) + up
   v1 = (v - vp) + vp;
   u2 = u - u1
   v2 = v - v1
   rh = u*v
   rl = (((u1 * v1 - rh) + u1 * v2) + u2 * v1) + u2 * v2
   return rh, rl 

def check_dekker(u,v):
   h, l = dekker (u,v)
   U = u.exact_rational()
   V = v.exact_rational()    
   H = h.exact_rational()
   L = l.exact_rational()
   return H+L == U*V

def check_dekker_all(K=10^6,rnd='RNDN'):
   R = RealField(64,rnd=rnd)
   for k in range(K):
      u = R.random_element()
      v = R.random_element()
      if not check_dekker(u,v):
         print (get_hex(u), get_hex(v))
         raise ValueError

from functools import cmp_to_key

def cmp(x,y):
   if x[2] < y[2]:
      return int(-1)
   if x[2] > y[2]:
      return int(1)
   # now x[2] = y[2]
   if x[1] < y[1]:
      return int(-1)
   if x[1] > y[1]:
      return int(1)
   if x[0] < y[0]:
      return int(-1)
   if x[0] > y[0]:
      return int(1)
   return int(0)

def cmpneg(x,y):
   if x[2] > y[2]:
      return int(-1)
   if x[2] < y[2]:
      return int(1)
   # now x[2] = y[2]
   if x[0] < y[0]:
      return int(-1)
   if x[0] > y[0]:
      return int(1)
   if x[1] < y[1]:
      return int(-1)
   if x[1] > y[1]:
      return int(1)
   return int(0)

# statall("/tmp/log")
# ([((9223372036854775808, -64), (18446744073709551616, 0))],
# [((-18446744073709551615, 0), (-13306513097844322491, -64))])
def statall(f):
   f = open(f,"r")
   l = []
   while true:
      s = f.readline()
      if s=='':
         break
      s = s.split(" ")
      assert len(s) == 5, "len(s) == 5"
      t0 = ZZ(s[0])
      t1 = ZZ(s[1])
      e = ZZ(s[2])
      assert ZZ(s[3])==64, "s[3]==64"
      assert ZZ(s[4])==64, "s[4]==64"
      l.append((t0,t1,e))
   f.close()
   lpos = [x for x in l if x[0]>0]
   lneg = [x for x in l if x[0]<0]
   assert len(l) == len(lpos) + len(lneg)
   lpos.sort(key=cmp_to_key(cmp))
   lpos2 = []
   for t0,t1,e in lpos:
      if lpos2==[]:
         lpos2 = [((t0,e),(t1,e))]
      else:
         t1old,e1old = lpos2[-1][1]
         if t1old*2^e1old > t0*2^e:
            print ((t1old,e1old), (t0, e))
         assert t1old*2^e1old <= t0*2^e, "t1old*2^e1old <= t0*2^e"
         if t1old*2^e1old == t0*2^e:
            lpos2[-1] = (lpos2[-1][0],(t1,e))
         else:
            lpos2.append(((t0,e),(t1,e)))
   lpos = lpos2
   lneg.sort(key=cmp_to_key(cmpneg))
   lneg2 = []
   for t0,t1,e in lneg:
      if lneg2==[]:
         lneg2 = [((t0,e),(t1,e))]
      else:
         t1old,e1old = lneg2[-1][1]
         if (t1old-1)*2^e1old > (t0-1)*2^e:
            print ((t1old,e1old), (t0, e))
         assert (t1old-1)*2^e1old <= (t0-1)*2^e, "(t1old-1)*2^e1old <= (t0-1)*2^e"
         if (t1old-1)*2^e1old == (t0-1)*2^e:
            lneg2[-1] = (lneg2[-1][0],(t1,e))
         else:
            lneg2.append(((t0,e),(t1,e)))
   lneg = lneg2
   return lpos,lneg

# given l a list of worst cases in [0,1), extend it to worst cases in [1,16384)
def extend_wc(l):
   extra = []
   R = RealField(64)
   for k in range(1,16384):
      lnew = []
      for x in l:
         assert x > 0, "x > 0"
         y = R(k)+x
         if y.exact_rational() == k + x.exact_rational():
            lnew.append(y)
            extra.append(y)
      l = lnew
   return extra  

# given l a list of worst cases in (-1,0], extend it to worst cases in (-16384,-1]
def extend_wc_neg(l):
   extra = []
   R = RealField(64)
   for k in range(-1,-16446,-1):
      lnew = []
      for x in l:
         assert x < 0, "x < 0"
         y = R(k)+x
         if y.exact_rational() == k + x.exact_rational():
            lnew.append(y)
            extra.append(y)
      l = lnew
   return extra  

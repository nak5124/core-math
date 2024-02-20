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
   print ("static const long double T2fast[32][2] = {")
   R = RealField(64)
   R32 = RealField(32)
   maxerr = 0
   for i in range(2^5):
      x = n(2^(i/2^5), 200)
      h = R(R32(x))
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
   print ("static const long double T1fast[32][2] = {")
   R = RealField(64)
   R32 = RealField(32)
   maxerr = 0
   for i in range(2^5):
      x = n(2^(i/2^10), 200)
      h = R(R32(x))
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
   print ("static const long double T0fast[32][2] = {")
   R = RealField(64)
   R32 = RealField(32)
   maxerr = 0
   for i in range(2^5):
      x = n(2^(i/2^15), 200)
      h = R(R32(x))
      l = R(x-h.exact_rational())
      print ("   {" + get_hex(h) + "L, " + get_hex(l) + "L},")
      maxerr = max(maxerr,abs(x-h.exact_rational()-l.exact_rational()))
   print ("};")
   print (log(maxerr)/log(2.))

# return the 'ulp' of the interval x, i.e., max(ulp(t)) for t in x
# this internal routine is used below
def RIFulp(x):
   return max(x.lower().ulp(),x.upper().ulp())

# analyze_P()
# err1= -115.999994496566
# err2= -97.9999944965657
# err3= -79.9999944965657
# err4= -80.0000000000000
# err5= -126.999984741211
# rel. err= -78.9472649184332
# max l= 1.08420217248550444e-19
def analyze_P():
   err0 = 2^-83.748 # absolute error
   R = RealField(64)
   RI = RealIntervalField(64)
   x = RI(-2^-16,2^-16)
   p4 = RI(R("0x1.3b2ab70cf131bd7ep-7",16))
   p3 = RI(R("0x1.c6b08d6835c26dep-5",16))
   p2 = RI(R("0x1.ebfbdff82c58ea86p-3",16))
   p1 = RI(R("0x1.62e42fefa39ef358p-1",16))
   p0 = RI(1)
   # y = p[4] * x + p[3]
   y = p4*x+p3
   err1 = (RIFulp(p4*x)+RIFulp(y))*x.abs().upper()^3
   print ("err1=", log(err1)/log(2.))
   # y = y * x + p[2]
   yin = y
   y = y*x+p2
   err2 = (RIFulp(yin*x)+RIFulp(y))*x.abs().upper()^2
   print ("err2=", log(err2)/log(2.))
   # *h = y * x + p[1]
   h = y*x+p1
   err3 = (RIFulp(y*x)+RIFulp(h))*x.abs().upper()
   print ("err3=", log(err3)/log(2.))
   # h = *h * x
   h = h*x
   err4 = RIFulp(h)
   print ("err4=", log(err4)/log(2.))
   # fast_two_sum (h, l, p[0], *h)
   h = p0+h
   u = RIFulp(h)
   l = RI(-u,u)
   err5 = 2^-127*h.abs().upper()
   print ("err5=", log(err5)/log(2.))
   # convert err0 into relative error
   err0 = err0*(h+l).abs().upper()
   err = err0+err1+err2+err3+err4+err5
   # convert into relative error
   err = err/(h+l).abs().lower()
   print ("rel. err=", log(err)/log(2.))
   print ("max l=", l.abs().upper())
   
# split binade [2^(e-1),2^e) into k blocks
def doit_bacsel(e,k):
   t0 = 2^63
   t1 = 2^64
   h = ceil((t1-t0)/k)
   for i in range(k):
      u0 = t0+h*i
      u1 = min(t0+h*(i+1),t1)
      print ("./doit.sh " + str(u0) + " " + str(u1) + " 64 " + str(e) + " 64 20 64")

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
         if t1old*2^e1old > t0*2^e:
            print ((t1old,e1old), (t0, e))
         assert (t1old-1)*2^e1old <= (t0-1)*2^e, "(t1old-1)*2^e1old <= (t0-1)*2^e"
         if (t1old-1)*2^e1old == (t0-1)*2^e:
            lneg2[-1] = (lneg2[-1][0],(t1,e))
         else:
            lneg2.append(((t0,e),(t1,e)))
   lneg = lneg2
   return lpos,lneg

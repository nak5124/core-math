# generate bacsel input for 2^(e-1) <= x < 2^e
def doit(e):
   if e <= -1022:
      n = 1074+e
   else:
      n = 53
   t0 = 2^(n-1)
   t1 = 2^n
   ein = e
   nthreads = 64
   t=20
   s = "./doit.sh " + str(t0) + " " + str(t1) + " " + str(n) + " " + str(ein)
   s += " " + str(nthreads) + " " + str(t)
   print (s)

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

# l=statall("out")
# [((1, -1074), (4503599627370496, -1074)), ((4503599627370496, -1026), (9007199254740992, -1023)), ((4503599627370496, -1021), (9007199254740992, -1019)), ((4503599627370496, -1017), (9007199254740992, -1013)), ((4503599627370496, -1011), (9007199254740992, -1010)), ((4503599627370496, -1007), (9007199254740992, -1006)), ((4503599627370496, -1004), (9007199254740992, -999)), ((4503599627370496, -997), (9007199254740992, -997)), ((4503599627370496, -995), (9007199254740992, -995)), ((4503599627370496, -993), (9007199254740992, -989)), ((4503599627370496, -106), (9007199254740992, -53))]
# the last entry ((4503599627370496, -106), (9007199254740992, -53)) means
# that all works cases were checked from 4503599627370496/2^106 ~ 5.5e-17
# to 9007199254740992/2^53 = 1
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

# generate hard-to-round cases in [1,2^52) from those in (0,1)
def gen_all(l):
   l = [x for x in l if x==x and 0 < x < 1]
   L = []
   for e in [1..53]: # hard-to-round cases in [2^(e-1),2^e)
      # ulp(x) = 2^(e-53) for x in [2^(e-1),2^e)
      l = [x for x in l if x.exact_rational()/2^(e-53) in ZZ]
      print ("e=", e, "remains", len(l))
      if len(l)==1: break # only remains 0.5
      x0 = 2^(e-1)
      while x0+1<2^e:
         L = L + [x0+x for x in l]
         x0 += 1
   return L


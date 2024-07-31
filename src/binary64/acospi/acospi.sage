#load("../../generic/support/common.sage")
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
# [((-9007199254740991, -53), (-7074237752028440, -105)), ((7074237752028441, -106), (9007199254740992, -53))]
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
      assert t0.nbits() == n, "t0.nbits() == n"
      assert (t1-1).nbits() == n, "(t1-1).nbits() == n"
      l.append((t0,t1,e-n))
   f.close()
   l.sort(key=cmp_to_key(cmp))
   l2 = []
   for t0,t1,e in l:
      if l2==[]:
         l2 = [((t0,e),(t1,e))]
      elif t0 > 0:
         t1old,e1old = l2[-1][1]
         if t1old*2^e1old > t0*2^e:
            print ((t1old,e1old), (t0, e))
         assert t1old*2^e1old <= t0*2^e, "t1old*2^e1old <= t0*2^e"
         if t1old*2^e1old == t0*2^e:
            l2[-1] = (l2[-1][0],(t1,e))
         else:
            l2.append(((t0,e),(t1,e)))
      else: # t0 < 0
         t1old,e1old = l2[-1][1]
         if (t1old-1).nbits() == t1old.nbits():
            next_x = t1old*2^e1old
         else:
            next_x = (2*t1old-1)*2^(e1old-1)
         if next_x > t0*2^e:
            print ((t1old,e1old), (t0, e))
         assert next_x <= t0*2^e, "next_x <= t0*2^e"
         if next_x == t0*2^e:
            l2[-1] = (l2[-1][0],(t1,e))
         else:
            l2.append(((t0,e),(t1,e)))
   l = l2
   return l

# from a list l of exceptional cases, print two tables:
# one with (x,h) where h is the rounding to nearest of f(x)
# one with a char l (-1,0,1) giving the direction of rounding
def format_exceptions(l):
   print ("static const double exceptions[EXCEPTIONS][2] = {")
   for x in l:
      X = x.exact_rational()
      h = RR(n(acospi(X),200))
      print ("    {" + get_hex(x) + ", " + get_hex(h) + "},")
   print ("  };")
   print ("static const char exceptions_rnd[EXCEPTIONS] = {")
   for x in l:
      X = x.exact_rational()
      h = RR(n(acospi(X),200))
      l = RR(n(acospi(X)-h.exact_rational(),200))
      if l==0:
         l=0
      elif l>0:
         l=1
      else:
         l=-1
      print ("    " + str(l) + ", /* " + get_hex(x) + " */")
   print ("  };")

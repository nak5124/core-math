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
# [((8109328451921441, -81), (9007199254740992, -53))]
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

# given l a list of worst-cases for 0 < x < 1, extend them
# for 1 < x < 2^52 (except exact cases)
def extend_wc(l):
   ll = []
   l = [x.exact_rational() for x in l]
   for e in [1..52]:
      # extend in 2^(e-1) <= x < 2^e
      for X in l:
         x = 2^(e-1)+X
         if 2*x in ZZ:
            continue
         if RR(x).exact_rational() == x:
            # for j in range(2^(e-1),2^e):
            for j in range(2^(e-1),2^e):
               ll.append(RR(j+X))
   return ll
            
# generate hard-to-round cases in [1,2^52) from those in (0,1)
def gen_all(l):
   l = [x for x in l if x==x and 0 <= x < 1 and x.sign_mantissa_exponent()[0]==1]
   l = [x for x in l if identical_bits(cospi,x)>=42] + [RR(0)]
   L = []
   for e in [1..53]: # hard-to-round cases in [2^(e-1),2^e)
      # ulp(x) = 2^(e-53) for x in [2^(e-1),2^e)
      l = [x for x in l if x.exact_rational()/2^(e-53) in ZZ]
      print ("e=", e, "remains", len(l))
      if len(l)==2: break # 0 and 0.5
      x0 = 2^(e-1)
      while x0+1<2^e:
         L = L + [x0+x for x in l]
         x0 += 1
   return L

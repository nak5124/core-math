#load("../../generic/support/common.sage")
def gen_tabi():
   print ("static const double tab_i[128][2] = {")
   for i in range(-127,128):
      x = 2^(i/128)
      h = RR(n(x,100))
      H = h.exact_rational()
      l = RR(n(x-H,200))
      print ("{ " + get_hex(h) + ", " + get_hex(l) + " }, /* " + str(i) + " */")
   print ("};")

def ph(h):
   h = h.exact_rational()
   p0 = 1
   p1 = RR("0x1.62e42fefa39efp-1",16).exact_rational()
   p1 += RR("0x1.b2ca0bb577094p-56",16).exact_rational()
   p2 = RR("0x1.ebfbdff82c58fp-3",16).exact_rational()
   p3 = RR("0x1.c6b08d70484c1p-5",16).exact_rational()
   p4 = RR("0x1.3b2ab6fb663a2p-7",16).exact_rational()
   p5 = RR("0x1.5d881a764d899p-10",16).exact_rational()
   p6 = RR("0x1.430bba9c70dddp-13",16).exact_rational()
   yh = p6
   yh = p5 + yh * h
   yh = p4 + yh * h
   yh = p3 + yh * h
   yh = p2 + yh * h
   yh = p1 + yh * h
   return yh

def check_fast_path(x,e,yh,yl):
   X = x.exact_rational()
   YH = yh.exact_rational()
   YL = yl.exact_rational()
   return n(exp(X)-2^e*(YH+YL))

# h=RR("0x1.e1af735b65cp-8",16) # x=0x1.011dd7222b8d2p+1
def check_accurate_poly(h):
   H = h.exact_rational()
   p0 = 1
   p1h = RR("0x1.62e42fefa39efp-1",16)
   p1l = RR("0x1.abc9e3b397ebp-56",16)
   p1 = p1h.exact_rational() + p1l.exact_rational()
   p2h = RR("0x1.ebfbdff82c58fp-3",16)
   p2l = RR("-0x1.5e43a5429b326p-57",16)
   p2 = p2h.exact_rational() + p2l.exact_rational()
   p3h = RR("0x1.c6b08d704a0cp-5",16)
   p3l = RR("-0x1.d331600cee073p-59",16)
   p3 = p3h.exact_rational() + p3l.exact_rational()
   p4h = RR("0x1.3b2ab6fba4e77p-7",16)
   p4l = RR("0x1.4fb30e5c2c8bcp-62",16)
   p4 = p4h.exact_rational() + p4l.exact_rational()
   p5 = RR("0x1.5d87fe78a6731p-10",16).exact_rational()
   p6 = RR("0x1.430912f86bfb8p-13",16).exact_rational()
   p7 = RR("0x1.ffcbfc58b51c9p-17",16).exact_rational()
   p8 = RR("0x1.62c034be4ffd9p-20",16).exact_rational()
   p9 = RR("0x1.b523023e3d552p-24",16).exact_rational()
   y = p8+H*p9
   y = p7+H*y
   y = p6+H*y
   y = p5+H*y
   y = p4+H*y
   y = p3+H*y
   y = p2+H*y
   y = p1+H*y
   y = p0+H*y
   e = n(2^H-y,200)
   return e, log(abs(e))/log(2.)

def gen_bacsel_aux(X0,X1,e,sign):
   scale = 2^(53-e)
   t0 = ZZ(X0.exact_rational()*scale)
   t1 = ZZ(X1.exact_rational()*scale)
   if sign==1:
      print ("./doit2.sh " + str(t0) + " " + str(t1) + " " + str(e))
   else:
      print ("./doit2.sh " + str(-t1+1) + " " + str(-t0+1) + " " + str(e))
   
def gen_bacsel(emin=0,emax=10,group=10,sign=1):
   assert sign in [1,-1]
   R = RealField(53,rnd='RNDZ')
   for e in [emin..emax]:
      x0 = R(2^(e-1))
      x1 = R(2^e)
      x10 = x1
      n = 0
      X0 = x0
      while true:
         x00 = x0
         if exp(sign*x0)>=2^1024.:
            break
         if exp(sign*x0)<2^-1075.:
            break
         while exp(x0).ulp() != exp(x1).ulp():
            x2 = (x0+x1)/2
            if x2 in [x0,x1]:
               break
            if exp(x2).ulp()==exp(x0).ulp():
               x0 = x2
            else:
               x1 = x2
         x0 = x00
         X1 = x1
         n += 1
         if n==group:
            gen_bacsel_aux(X0,X1,e,sign)
            X0 = X1
            n = 0
         if x1==x10:
            break
         x0 = x1
         x1 = x10
      if n>0:
         gen_bacsel_aux(X0,X1,e,sign)

def get_hl(x):
   X = x.exact_rational()
   h = RR(n(exp(X),200))
   H = h.exact_rational()
   l = RR(n(exp(X)-H,200))
   L = l.exact_rational()
   if abs(l)==h.ulp()/2 or l==0:
      if n(exp(X)-H-L,200)>0:
         l = l.nextabove()
      else:
         l = l.nextbelow()
   return h, l

def gen_exceptions_aux(x):
   h, l = get_hl (x)
   print ("      if (x == " + get_hex(x) + ")")
   if l > 0:
      print ("        return " + get_hex(h) + " + " + get_hex(l) + ";")
   else:
      s = get_hex(l)
      print ("        return " + get_hex(h) + " - " + s[1:] + ";")

# generate the code for exceptional values in l
# l=[RR("0x1.62f71c4656b61p-1",16),RR("0x1.333a83013057ep+2",16),RR("0x1.83d4bcdebb3f4p+2",16)]
# gen_exceptions(l)
def gen_exceptions(l):
   # first sort l by value of e
   T = dict()
   for x in l:
      n = asuint64(x)
      e = n % 64
      if not e in T.keys():
         T[e] = []
      T[e].append(x)
   # now generate code
   print ("  switch (e) {")
   for e in T.keys():
      le = T[e]
      print ("    case "+str(e)+":")
      for x in le:
         gen_exceptions_aux(x)
      print ("      break;")
   print ("  };")
   
# output parameters for bacsel for subnormal range 2^(e-1) <= exp(x) < 2^e
def check_subnormal(emin=-1073,emax=-1022,nthreads=64):
   for e in [emin..emax]:
      assert -1073 <= e <= -1022
      # 2^(e-1) <= exp(x) thus log(2^(e-1)) <= x
      R = RealField(53,rnd='RNDU')
      xmin = log(R(2^(e-1)))
      xmax = log(R(2^(e)))
      Xmin = xmin.exact_rational()
      Xmax = xmax.exact_rational()
      emin = 1+floor(log(abs(Xmin))/log(2))
      emax = 1+floor(log(abs(Xmin))/log(2))
      assert emin==emax, "emin==emax"
      t0 = 2^53*Xmin/2^emin
      assert t0 in ZZ, "t0 in ZZ"
      t0 = ZZ(t0)
      t1 = 2^53*Xmax/2^emax
      assert t1 in ZZ, "t1 in ZZ"
      t1 = ZZ(t1)
      nn = 1074+e
      # around t0,t1 we use d=0 to avoid failures in bacsel
      margin = 2^20
      t0a = min(t0+margin,t1)
      t1a = max(t0,t1-margin)
      if t0<t0a:
         print ("./bacsel -rnd_mode all -prec 128 -n 53 -nn " + str(nn) + " -m 44 -t 20 -t0 " + str(t0) + " -t1 " + str(t0a) + " -d 0 -alpha 1 -v -e_in " + str(emin) + " -nthreads " + str(nthreads))
      if t0a<t1a:
         print ("./bacsel -rnd_mode all -prec 128 -n 53 -nn " + str(nn) + " -m 44 -t 20 -t0 " + str(t0a) + " -t1 " + str(t1a) + " -d 2 -alpha 2 -v -e_in " + str(emin) + " -nthreads " + str(nthreads))
      if t1a<t1:
         print ("./bacsel -rnd_mode all -prec 128 -n 53 -nn " + str(nn) + " -m 44 -t 20 -t0 " + str(t1a) + " -t1 " + str(t1) + " -d 0 -alpha 1 -v -e_in " + str(emin) + " -nthreads " + str(nthreads))

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

from functools import cmp_to_key

# lpos,lneg = statall("/tmp/log")
# lpos = [((4503599627370496, -52), (6243314768165360, 10))]
# lneg = [((-6554261109157969, 10), (-4503599627370495, -54))]
def statall(f):
   f = open(f,"r")
   l = []
   while true:
      s = f.readline()
      if s=='':
         break
      s = s.split(" ")
      assert len(s) == 3
      t0 = ZZ(s[0])
      t1 = ZZ(s[1])
      e = ZZ(s[2])
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

def statall_denorm(f):
   f = open(f,"r")
   l = []
   while true:
      s = f.readline()
      if s=='':
         break
      s = s.split(" ")
      assert len(s) == 4
      t0 = ZZ(s[0])
      t1 = ZZ(s[1])
      e = ZZ(s[2])
      nn = ZZ(s[3])
      l.append((t0,t1,e,nn))
   f.close()
   # check each range
   for t0,t1,e,nn in l:
      # check 2^(-1074+nn-1) <= exp(t0*2^(e-53))
      assert 2^(-1074+nn-1) <= exp(t0*2^(e-53)), "assert 1"
      # check exp((t0-1)*2^(e-53)) < 2^(-1074+nn-1)
      assert exp((t0-1)*2^(e-53)) < 2^(-1074+nn-1), "assert 2"
      # check 2^(-1074+nn) <= exp(t1*2^(e-53))
      assert 2^(-1074+nn) <= exp(t1*2^(e-53)), "assert 3"
      # check exp((t1-1)*2^(e-53)) < 2^(-1074+nn)
      assert exp((t1-1)*2^(e-53)) < 2^(-1074+nn), "assert 4"
   # now check there is no hole
   l.sort(key=cmp_to_key(cmp))
   for i in range(len(l)-1):
      assert l[i][3] == i+1, "l[i][3] == i+1"
      assert l[i+1][3] == i+2, "l[i+1][3] == i+2"
      assert l[i][2] == 10, "l[i][2] == 10"
      assert l[i+1][2] == 10, "l[i+1][2] == 10"
      assert l[i][1] == l[i+1][0], "l[i][1] == l[i+1][0]"
   return [(l[0][0],l[-1][1],10)]
         
def muldd(xh,xl,ch,cl):
   ahhh = ch*xh
   l = (ch*xl + cl*xh) + fma(ch, xh, -ahhh)
   return ahhh, l

def fasttwosum(x,y):
   s = x+y
   z = s-x
   return s, y-z

def fastsum(xh,xl,yh,yl):
   sh, sl = fasttwosum(xh, yh)
   return sh, (xl + yl) + sl

def FastFMA_DW(ah,al,bh,bl,ch,cl):
   dh = fma (ah, bh, ch)
   t = ch - dh
   e = fma (ah, bh, t)
   f = e + cl
   g = ah * bl + f
   l = al * bh + g
   return dh, l

def opolydd(xh,xl,n,c):
   h, l = c[n-1]
   for i in range(n-2,-1,-1):
      h, l = FastFMA_DW(xh,xl,h,l,c[i][0],c[i][1])
   return h, l

def check_accurate(x):
   if type(x)==type(''):
      x = RR(x,16)
   global T0, T1
   for r in 'NZUD':
      R = RealField(53,rnd='RND'+r)
      x = R(x)
      s = R("0x1.71547652b82fep+12",16)
      jt = round(x*s)
      t = R(jt)
      i0 = (jt>>6)&0x3f
      i1 = jt&0x3f
      ie = jt>>12
      t0h = R(T0[i0][1],16)
      t0l = R(T0[i0][0],16)
      t1h = R(T1[i1][1],16)
      t1l = R(T1[i1][0],16)
      th, tl = muldd(t0h,t0l, t1h,t1l)
      l2h = R("0x1.62e42ffp-13",16)
      l2l = R("0x1.718432a1b0e26p-47",16)
      l2ll = R("0x1.9ff0342542fc3p-102",16)
      dx = x - l2h*t
      dxl = l2l*t
      dxll = l2ll*t + fma(l2l,t,-dxl)
      dxh, dxl = fasttwosum (dx, dxl)
      dxl += dxll
      ch = [("0x1p+0","0"),("0x1p-1","0x1.712f72ecec2cfp-99"),("0x1.5555555555555p-3","0x1.5555555554d07p-57"),("0x1.5555555555555p-5","0x1.55194d28275dap-59"),("0x1.1111111111111p-7","0x1.12faa0e1c0f7bp-63"),("0x1.6c16c16da6973p-10","-0x1.4ba45ab25d2a3p-64"),("0x1.a01a019eb7f31p-13","-0x1.9091d845ecd36p-67")]
      ch = [(R(x[0],16),R(x[1],16)) for x in ch]
      fh, fl = opolydd(dxh,dxl,7,ch)
      fh,fl = muldd(dxh,dxl,fh,fl)
      fh,fl = muldd(fh,fl,th,tl)
      fh,fl = fastsum(th,tl,fh,fl)
      fh,fl = fasttwosum(fh,fl)
      # print ("fh=", get_hex(fh))
      z = fh*2^ie
      # print (get_hex(z), get_hex(exp(x)))
      if z!=exp(x): return false
   return true

def get_i0i1(x):
   R = x.parent()
   s = R("0x1.71547652b82fep+12",16)
   jt = round(x*s)
   i0 = (jt>>6)&0x3f
   i1 = jt&0x3f
   return i0,i1

# width=1: remains 28
# width=2: remains 20
# width=2: remains 17
# width=4: remains 15
# width=5: remains 14
# optimize the number of exceptional cases in the accurate path
def optimize(width=1):
   global T0, T1, Exceptions
   exceptions = [RR(x,16) for x in Exceptions]
   L = [x for x in exceptions if check_accurate(x)==false]
   n = len(L)
   print ("n=", n)
   Changed = true
   W = []
   for w in [1..width]:
      W.append(w)
      W.append(-w)
   oldn = n # current number of exceptions
   while Changed:
      Changed = false
      for i0 in range(2^6):
         li0 = [x for x in L if get_i0i1(x)[0]==i0]
         if li0!=[]:
            # try to modify T0[i0]
            l,h = T0[i0]
            y = RR(l,16)
            changed = false
            for d in W:
               T0[i0] = (get_hex(y+d*y.ulp()),h)
               li0 = [x for x in exceptions if check_accurate(x)==false]
               n = len(li0)
               if n<oldn:
                  print ("i0=", i0, "reduced from", oldn, "to", n)
                  oldn = n
                  L = li0
                  changed = Changed = true
                  break
            if changed==false:
               T0[i0] = l,h # restore initial values
      for i1 in range(2^6):
         li1 = [x for x in L if get_i0i1(x)[1]==i1]
         if li1!=[]:
            # try to modify T1[i1]
            l,h = T1[i1]
            y = RR(l,16)
            changed = false
            for d in W:
               T1[i1] = (get_hex(y+d*y.ulp()),h)
               li1 = [x for x in exceptions if check_accurate(x)==false]
               n = len(li1)
               if n<oldn:
                  print ("i1=", i1, "reduced from", oldn, "to", n)
                  oldn = n
                  L = li1
                  changed = Changed = true
                  break
            if changed==false:
               T1[i1] = l,h # restore initial values

# print differences between original and new versions
# %runfile exp.sage
# optimize(width=5)
# newT0,newT1=T0,T1
# %runfile exp.sage # reload the original tables
# print_diff(newT0,newT1)
def print_diff(newT0,newT1):
   global T0,T1
   for i0 in range(2^6):
      if T0[i0]!=newT0[i0]:
         assert T0[i0][1]==newT0[i0][1]
         print ("i0=", i0, T0[i0][0], "changed to", newT0[i0][0])
   for i1 in range(2^6):
      if T1[i1]!=newT1[i1]:
         assert T1[i1][1]==newT1[i1][1]
         print ("i1=", i1, T1[i1][0], "changed to", newT1[i1][0])

Exceptions = ["0x1p-53", "0x1.ffff7fffe0001p-36", "0x1.8cd99fffb319cp-34", "0x1.413f33fe6ce07p-31", "0x1.f8e165f8388f7p-30", "0x1.282aae3169ee7p-20", "0x1.e5ba92aa9b1efp-20", "0x1.ba07d73250de7p-14", "0x1.6a4d1af9cc989p-8", "0x1.5a75293a5dcdap-6", "0x1.42ea46949b3c7p-5", "0x1.7c8bb0cf5d16p-5", "0x1.0948d39a41695p-3", "0x1.a065fefae814fp-3", "0x1.f6e4c3ced7c72p-3", "0x1.1a0408712e00ap-2", "0x1.bcab27d05abdep-2", "0x1.005ae04256babp-1", "0x1.273c188aa7b14p+2", "0x1.83d4bcdebb3f4p+2", "0x1.08f51434652c3p+4", "0x1.1d5c2daebe367p+4", "0x1.c44ce0d716a1ap+4", "0x1.e07e71bfcf06fp+5", "0x1.f7216c4b435c9p+5", "0x1.54cd1fea7663ap+7", "0x1.d6479eba7c971p+8", "-0x1.4bd46601ae1efp-31", "-0x1.85e60704a3a9cp-30", "-0x1.ff171507f8ba5p-30", "-0x1.0401ae48409b5p-28", "-0x1.664716b68a409p-14", "-0x1.a2fefefd580dfp-13", "-0x1.ce3f638d0c742p-12", "-0x1.ceff32831e2c2p-12", "-0x1.33accae78b371p-11", "-0x1.d792b60084f92p-11", "-0x1.7fb235d76cce7p-8", "-0x1.1ff9b8e8b38bep-7", "-0x1.54511e930898cp-7", "-0x1.5c5ed0ec83666p-6", "-0x1.8c56ff5326197p-6", "-0x1.a4187f2ca71f9p-6", "-0x1.a8f783d749a8fp-4", "-0x1.bd44fdaed819fp-4", "-0x1.daf693d64fadap-4", "-0x1.290ea09e36479p-3", "-0x1.8aeb636f3ce35p-3", "-0x1.d3f3799439415p-3", "-0x1.ea16274b0109bp-3", "-0x1.22e24fa3d5cf9p-1", "-0x1.85068c07fbbf6p-1", "-0x1.bdc7955d1482cp-1", "-0x1.2a9cad9998262p+0", "-0x1.cc37ef7de7501p+0", "-0x1.02393d5976769p+1", "-0x1.65061daf79a78p+1", "-0x1.e8bdbfcd9144ep+3", "-0x1.8f80e06f3a04cp+4", "-0x1.59f038076039cp+6", "-0x1.981587ad4542fp+7","-0x1.07e05e6e03b83p-2","-0x1.1b3ad95690545p-13","-0x1.36ca6d2bbd54ep+6","-0x1.3f53708c0d3dp-13","-0x1.40a990279c2cp-6","-0x1.66cdaa8c7998dp-11","-0x1.8c100c37daa9ep-7","-0x1.daf3afee51991p+2","0x1.2f3832cad3d5fp+1","0x1.7a60ee15e3e9dp+6","0x1.aca7ae8da5a7bp+0","0x1.d086543694c5ap-1","0x1.293c6f3f245f9p-2"]

T0 = [
  ("0x0p+0", "0x1p+0"), ("-0x1.19083535b085cp-56", "0x1.02c9a3e778061p+0"),
  ("0x1.d73e2a475b466p-55", "0x1.059b0d3158574p+0"), ("0x1.186be4bb284fdp-57", "0x1.0874518759bc8p+0"),
  ("0x1.8a62e4adc610bp-54", "0x1.0b5586cf9890fp+0"), ("0x1.03a1727c57b55p-59", "0x1.0e3ec32d3d1a2p+0"),
  ("-0x1.6c51039449b3ap-54", "0x1.11301d0125b51p+0"), ("-0x1.32fbf9af1369ep-54", "0x1.1429aaea92dep+0"),
  ("-0x1.19041b9d78a76p-55", "0x1.172b83c7d517bp+0"), ("0x1.e5b4c7b4968e4p-55", "0x1.1a35beb6fcb75p+0"),
  ("0x1.e016e00a2643cp-54", "0x1.1d4873168b9aap+0"), ("0x1.dc775814a8496p-55", "0x1.2063b88628cd6p+0"),
  ("0x1.9b07eb6c70572p-54", "0x1.2387a6e756238p+0"), ("0x1.2bd339940e9dap-55", "0x1.26b4565e27cddp+0"),
  ("0x1.612e8afad1255p-55", "0x1.29e9df51fdee1p+0"), ("0x1.0024754db41d4p-54", "0x1.2d285a6e4030bp+0"),
  ("0x1.6f46ad23182e4p-55", "0x1.306fe0a31b715p+0"), ("0x1.32721843659a6p-54", "0x1.33c08b26416ffp+0"),
  ("-0x1.63aeabf42eaep-54", "0x1.371a7373aa9cbp+0"), ("-0x1.5e436d661f5dbp-56", "0x1.3a7db34e59ff7p+0"),
  ("0x1.ada0911f09ebap-55", "0x1.3dea64c123422p+0"), ("-0x1.ef3691c309278p-58", "0x1.4160a21f72e2ap+0"),
  ("0x1.89b7a04ef80dp-59", "0x1.44e086061892dp+0"), ("0x1.3c1a3b69062fp-56", "0x1.486a2b5c13cdp+0"),
  ("0x1.d4397afec42e2p-56", "0x1.4bfdad5362a27p+0"), ("-0x1.4b309d25957e3p-54", "0x1.4f9b2769d2ca7p+0"),
  ("-0x1.07abe1db13cabp-55", "0x1.5342b569d4f82p+0"), ("0x1.9bb2c011d93acp-54", "0x1.56f4736b527dap+0"),
  ("0x1.6324c054647acp-54", "0x1.5ab07dd485429p+0"), ("0x1.ba6f93080e65ep-54", "0x1.5e76f15ad2148p+0"),
  ("-0x1.383c17e40b497p-54", "0x1.6247eb03a5585p+0"), ("-0x1.bb60987591c34p-54", "0x1.6623882552225p+0"),
  ("-0x1.bdd3413b26456p-54", "0x1.6a09e667f3bcdp+0"), ("-0x1.bbe3a683c88aap-57", "0x1.6dfb23c651a2fp+0"),
  ("-0x1.16e4786887a9ap-55", "0x1.71f75e8ec5f74p+0"), ("-0x1.0245957316dd4p-54", "0x1.75feb564267c9p+0"),
  ("-0x1.41577ee04993p-55", "0x1.7a11473eb0187p+0"), ("0x1.05d02ba15797ep-56", "0x1.7e2f336cf4e62p+0"),
  ("-0x1.d4c1dd41532d8p-54", "0x1.82589994cce13p+0"), ("-0x1.fc6f89bd4f6bap-54", "0x1.868d99b4492edp+0"),
  ("0x1.6e9f156864b29p-54", "0x1.8ace5422aa0dbp+0"), ("0x1.5cc13a2e3976fp-55", "0x1.8f1ae99157736p+0"),
  ("-0x1.75fc781b57ebcp-57", "0x1.93737b0cdc5e5p+0"), ("-0x1.d185b7c1b85dp-54", "0x1.97d829fde4e5p+0"),
  ("0x1.c7c46b071f2bep-56", "0x1.9c49182a3f09p+0"), ("-0x1.359495d1cd532p-54", "0x1.a0c667b5de565p+0"),
  ("-0x1.d2f6edb8d41e3p-54", "0x1.a5503b23e255dp+0"), ("0x1.0fac90ef7fd32p-54", "0x1.a9e6b5579fdbfp+0"),
  ("0x1.7a1cd345dcc82p-54", "0x1.ae89f995ad3adp+0"), ("-0x1.2805e3084d708p-57", "0x1.b33a2b84f15fbp+0"),
  ("-0x1.5584f7e54ac3bp-56", "0x1.b7f76f2fb5e47p+0"), ("0x1.23dd07a2d9e84p-55", "0x1.bcc1e904bc1d2p+0"),
  ("0x1.11065895048dep-55", "0x1.c199bdd85529cp+0"), ("0x1.2884dff483cadp-54", "0x1.c67f12e57d14bp+0"),
  ("0x1.503cbd1e949d9p-56", "0x1.cb720dcef9069p+0"), ("-0x1.cbc3743797a9cp-54", "0x1.d072d4a07897cp+0"),
  ("0x1.2ed02d75b3706p-55", "0x1.d5818dcfba487p+0"), ("0x1.c2300696db533p-54", "0x1.da9e603db3285p+0"),
  ("-0x1.1a5cd4f184b5cp-54", "0x1.dfc97337b9b5fp+0"), ("0x1.39e8980a9cc9p-55", "0x1.e502ee78b3ff6p+0"),
  ("-0x1.e9c23179c2894p-54", "0x1.ea4afa2a490dap+0"), ("0x1.dc7f486a4b6b1p-54", "0x1.efa1bee615a27p+0"),
  ("0x1.9d3e12dd8a18dp-54", "0x1.f50765b6e454p+0"), ("0x1.74853f3a5931fp-55", "0x1.fa7c1819e90d8p+0")
]

T1 = [
  ("0x0p+0", "0x1p+0"), ("0x1.ae8e38c59c72cp-54", "0x1.000b175effdc7p+0"),
  ("-0x1.7b5d0d58ea8f4p-58", "0x1.00162f3904052p+0"), ("0x1.4115cb6b16a8dp-54", "0x1.0021478e11ce6p+0"),
  ("-0x1.d7c96f201bb2ep-55", "0x1.002c605e2e8cfp+0"), ("0x1.84711d4c35eap-54", "0x1.003779a95f959p+0"),
  ("-0x1.0484245243778p-55", "0x1.0042936faa3d8p+0"), ("-0x1.4b237da2025fap-54", "0x1.004dadb113dap+0"),
  ("-0x1.5e00e62d6b30ep-56", "0x1.0058c86da1c0ap+0"), ("0x1.a1d6cedbb948p-54", "0x1.0063e3a559473p+0"),
  ("-0x1.4acf197a00141p-54", "0x1.006eff583fc3dp+0"), ("-0x1.eaf2ea42391a6p-57", "0x1.007a1b865a8cap+0"),
  ("0x1.da93f90835f7ap-56", "0x1.0085382faef83p+0"), ("-0x1.6a79084ab0937p-55", "0x1.00905554425d4p+0"),
  ("0x1.86364f8fbe8f8p-54", "0x1.009b72f41a12bp+0"), ("-0x1.82e8e14e3110bp-55", "0x1.00a6910f3b6fdp+0"),
  ("-0x1.4f6b2a7609f72p-55", "0x1.00b1afa5abcbfp+0"), ("-0x1.e1a258ea8f71ap-56", "0x1.00bcceb7707ecp+0"),
  ("0x1.4362ca5bc26f2p-56", "0x1.00c7ee448ee02p+0"), ("0x1.095a56c919d02p-54", "0x1.00d30e4d0c483p+0"),
  ("-0x1.406ac4e81a646p-57", "0x1.00de2ed0ee0f5p+0"), ("0x1.b5a6902767e08p-54", "0x1.00e94fd0398ep+0"),
  ("-0x1.91b206085932p-54", "0x1.00f4714af41d3p+0"), ("0x1.427068ab22306p-55", "0x1.00ff93412315cp+0"),
  ("0x1.c1d0660524e08p-54", "0x1.010ab5b2cbd11p+0"), ("-0x1.e7bdfb3204be8p-54", "0x1.0115d89ff3a8bp+0"),
  ("0x1.843aa8b9cbbc6p-55", "0x1.0120fc089ff63p+0"), ("-0x1.34104ee7edae8p-56", "0x1.012c1fecd613bp+0"),
  ("-0x1.2b6aeb6176892p-56", "0x1.0137444c9b5b5p+0"), ("0x1.a8cd33b8a1bb3p-56", "0x1.01426927f5278p+0"),
  ("0x1.2edc08e5da99ap-56", "0x1.014d8e7ee8d2fp+0"), ("0x1.57ba2dc7e0c72p-55", "0x1.0158b4517bb88p+0"),
  ("0x1.b61299ab8cdb8p-54", "0x1.0163da9fb3335p+0"), ("-0x1.90565902c5f43p-54", "0x1.016f0169949edp+0"),
  ("0x1.70fc41c5c2d54p-55", "0x1.017a28af25567p+0"), ("0x1.4b9a6e145d76cp-54", "0x1.018550706ab62p+0"),
  ("-0x1.008eff5142bfap-56", "0x1.019078ad6a19fp+0"), ("-0x1.77669f033c7dep-54", "0x1.019ba16628de2p+0"),
  ("-0x1.09bb78eeead0ap-54", "0x1.01a6ca9aac5f3p+0"), ("0x1.371231477ece4p-54", "0x1.01b1f44af9f9ep+0"),
  ("0x1.5e7626621eb5ap-56", "0x1.01bd1e77170b4p+0"), ("-0x1.bc72b100828a4p-54", "0x1.01c8491f08f08p+0"),
  ("-0x1.ce39cbbab8bbep-57", "0x1.01d37442d507p+0"), ("0x1.16996709da2e5p-55", "0x1.01de9fe280ac8p+0"),
  ("-0x1.c11f5239bf536p-55", "0x1.01e9cbfe113efp+0"), ("0x1.e1d4eb5edc6b4p-55", "0x1.01f4f8958c1c6p+0"),
  ("-0x1.afb99946ee3fp-54", "0x1.020025a8f6a35p+0"), ("-0x1.8f06d8a148a32p-54", "0x1.020b533856324p+0"),
  ("-0x1.2bf310fc54eb8p-55", "0x1.02168143b0281p+0"), ("-0x1.c95a035eb4177p-54", "0x1.0221afcb09e3ep+0"),
  ("-0x1.491793e46834cp-54", "0x1.022cdece68c4fp+0"), ("-0x1.3e8d0d9c4908fp-56", "0x1.02380e4dd22adp+0"),
  ("-0x1.314aa16278aa4p-54", "0x1.02433e494b755p+0"), ("0x1.48daf888e965p-55", "0x1.024e6ec0da046p+0"),
  ("0x1.56dc8046821f4p-55", "0x1.02599fb483385p+0"), ("0x1.45b42356b9d46p-54", "0x1.0264d1244c719p+0"),
  ("-0x1.082ef51b61d7ep-56", "0x1.027003103b10ep+0"), ("0x1.2106ed0920a32p-56", "0x1.027b357854772p+0"),
  ("-0x1.fd4cf26ea5d11p-54", "0x1.0286685c9e059p+0"), ("-0x1.09f8775e78085p-54", "0x1.02919bbd1d1d8p+0"),
  ("0x1.64cbba902ca25p-58", "0x1.029ccf99d720ap+0"), ("0x1.4383ef231d208p-54", "0x1.02a803f2d170dp+0"),
  ("0x1.4a47a505b3a46p-54", "0x1.02b338c811703p+0"), ("0x1.e47120223467fp-54", "0x1.02be6e199c811p+0"),
]


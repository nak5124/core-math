def doit_bacsel_subnormal():
   for e in [-1073..-1022]:
      nn = 1074 + e # number of bits of output
      # deal with |exp2m1(x)| in [2^(e-1),2^e)
      # first deal with positive values
      a = RR(2^-1074)
      b = RR(1)
      while a.nextabove()!=b:
         c = (a+b)/2
         if n(c*log(2)) < 2^(e-1):
            a = c
         else:
            b = c
      x0 = b
      a = b
      b = RR(1)
      while a.nextabove()!=b:
         c = (a+b)/2
         if n(c*log(2)) < 2^e:
            a = c
         else:
            b = c
      x1 = b
      t0, e0, n0 = t_exp (x0)
      t1, e1, n1 = t_exp (x1)
      print_bacsel_pos (t0, e0, n0, t1, e1, n1, nn)

def t_exp(x):
   if abs(x)<2^-1022:
      m = ceil(x.exact_rational()*2^1074)
      n = m.nbits()
      assert n<=52, "n<=52"
      return m, -1074+n, n
   else:
      s,m,e = x.sign_mantissa_exponent()
      assert m.nbits()==53, "m.nbits()==53"
      return s*m, e+53, 53

def print_bacsel_pos (t0, e0, n0, t1, e1, n1, nn):
   assert n0<=n1, "n0<=n1"
   if n0<n1:
      print_bacsel_pos (t0, e0, n0, 2^n0, e0, n0, nn)
      print_bacsel_pos (2^n0, e0+1, n0+1, t1, e1, n1, nn)
   else:
      assert e0<=e1, "e0<=e1"
      if e0<e1:
         print_bacsel_pos (t0, e0, n0, 2^n0, e0, n0, nn)
         print_bacsel_pos (2^(n0-1), e0+1, n0, t1, e1, n1, nn)
      else:
         print ("./doit1.sh " + str(t0) + " " + str(t1) + " " + str(n0) + " " + str(e0) + " 64 10 " + str(nn) + " >> out")

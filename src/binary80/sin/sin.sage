# generate worst cases with sin(x) near +1/-1
def wc(out):
   f = open(out,"w")
   R = RealField(64)
   for e in [0..16384]:
      u = 2^(e-64)
      l = continued_fraction(u/(pi/2))
      bestq = 0
      for r in l.convergents():
         p,q = r.numer(),r.denom()
         # q*u ~ p*pi/2
         if q.nbits()>64:
            break
         if is_odd(p): # if p is odd, then sin(q*u) is near 1 or -1
            bestq = q
      if bestq != 0:
         x = R(bestq*u)
         f.write(get_hex(x)+"\n")
   f.close()
      



    

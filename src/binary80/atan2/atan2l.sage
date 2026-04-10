# print K hard-to-round cases in file f using the algorithm described in
# "The CORE-MATH Project",
# by Alexei Sibidanov, Paul Zimmermann, Stéphane Glondu
# Proceedings of the 29th IEEE Symposium on Computer Arithmetic (ARITH 2022)
# end of Section IIB
def worst_cases(f,K):
   f = open(f,"w")
   R = RealField(65)
   for k in range(K):
      z = R.random_element()
      Z = z.exact_rational()
      t = n(tan(Z),200)
      l = continued_fraction(t)
      for r in l.convergents():
         if r.numer().nbits()>64 or r.denom().nbits()>64:
            break
         oldr = r
      y = R(oldr.numer())
      x = R(oldr.denom())
      f.write(get_hex(y)+","+get_hex(x)+"\n")
   f.close()

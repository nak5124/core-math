# print K hard-to-round cases in file f using the algorithm described in
# "The CORE-MATH Project",
# by Alexei Sibidanov, Paul Zimmermann, Stéphane Glondu
# Proceedings of the 29th IEEE Symposium on Computer Arithmetic (ARITH 2022)
# end of Section IIB
# worst_cases("/tmp/out.wc",1000)
# 0x1.21e66e10d7631afp-4 0x1.3e0f2dea9e840618p-3 72
# worst_cases("/tmp/out.wc",1000,minm=64)
# 0x1.99a3f2d3e5372762p-1 0x1.bd46a1e2c476406cp-1 74
# worst_cases("/tmp/out.wc",1000,minm=70)
# -0x1.56fb80a640085438p+61 0x1.7c88fabc664ef938p+61 79
def worst_cases(f,K,wanted_m=40):
   maxm = 0
   f = open("/tmp/out"+str(wanted_m)+".wc","w")
   R = RealField(65)
   while K>0:
      z = R.random_element()
      Z = z.exact_rational()
      t = n(tan(Z),200)
      l = continued_fraction(t)
      for r in l.convergents():
         if r==0:
            continue
         if r.numer().nbits()>64 or r.denom().nbits()>64:
            break
         y = R(r.numer())/2^64
         x = R(r.denom())/2^64
         m = identical_bits_atan2(y,x)
         if m==wanted_m:
            f.write(get_hex(y)+","+get_hex(x)+"\n")
            K -= 1
            break
   f.close()

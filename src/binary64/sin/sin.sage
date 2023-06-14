# compute table for 1/(2pi)
def computeT():
   x = n(1/(2*pi),2000)
   m = ZZ(floor(x.exact_rational()*2^1216))
   l = m.digits(2^64)
   print ("static const uint64_t T[" + str(len(l)) + "] = {")
   for i in range(len(l)-1,-1,-1):
      print ("   " + hex(l[i]) + ",")
   print ("};")

# given x a 128-bit float, print it
def out_dint(x):
   s, m, e = x.sign_mantissa_exponent()
   while m!=0 and m.nbits()<128:
      m *= 2
      e -= 1
   h, l = divmod(m,2^64)
   print ("  {.hi = " + hex(h) + ", .lo = " + hex(l) + ", .ex = " + str(e+128) + ", .sgn=" + str((1-s)/2) + "},")

# compute table of sin(2*pi*i/2^8) for 0 <= i < 256
def computeS():
   R = RealField(128)
   print ("static const dint64_t S[256] = {")
   for i in range(256):
      s = n(sin(2*pi*i/2^8), 512)
      s = R(s)
      out_dint(s)
   print ("};")

# compute table of cos(2*pi*i/2^8) for 0 <= i < 256
def computeC():
   R = RealField(128)
   print ("static const dint64_t C[256] = {")
   for i in range(256):
      s = n(cos(2*pi*i/2^8), 512)
      s = R(s)
      out_dint (s)
   print ("};")

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

# compute table of sin(2*pi*i/2^11) for 0 <= i < 256
def computeS():
   R = RealField(128)
   print ("static const dint64_t S[256] = {")
   for i in range(256):
      s = n(sin(2*pi*i/2^11), 512)
      s = R(s)
      out_dint(s)
   print ("};")

# compute table of cos(2*pi*i/2^11) for 0 <= i < 256
def computeC():
   R = RealField(128)
   print ("static const dint64_t C[256] = {")
   for i in range(256):
      s = n(cos(2*pi*i/2^11), 512)
      s = R(s)
      out_dint (s)
   print ("};")

# implement Algorithm 2.1 from https://www.vinc17.net/research/papers/these.pdf
# finds the smallest distance between y=a*x-b and integers, 0 <= x < N
# return d the smallest distance, and x the corresponding value
def AlgoLefevre(a,b,N):
   x = frac(a)
   y = 1-frac(a)
   d = frac(b)
   u = v = 1
   j = 0
   # invariant: d = j*frac(a) - frac(b) mod 1
   while true:
      if d<x:
         while x<y:
            if u+v>=N:
               break
            y = y-x
            u = u+v
         if u+v>=N:
            break
         x = x-y
         v = v+u
      else:
         d = d-x
         j += v
         while y<x:
            if u+v>=N:
               break
            x = x-y
            v = v+u
         if u+v>=N:
            break
         y = y-x
         u = u+v
      if u+v>=N:
         break
   return d, j

# return x with smallest value of |x/(2pi) cmod 1| for 2^(e-1) <= x < 2^e
def search(e):
   R = RealField(e+200)
   N = 2^52
   b = R(2^(e-1)/(2*pi))
   b = frac(-b)
   if b < 0:
     b = 1+b
   assert b > 0, "b > 0"
   a = R(2^(e-53)/(2*pi))
   a = frac(a)
   assert a > 0, "a > 0"
   oldx = -1
   d, x = AlgoLefevre(a,b,N)
   d, x = AlgoLefevre(a,b+R(d),N)
   x = RR(2^(e-1)+x*2^(e-53))
   return x

def search_all():
   for e in range(1024,1,-1):
      x = search(e)
      print (get_hex(x) + ' # ' + str(e))

# for 2^(e-1) <= x < 2^e
# sin(x) is monotonous between (k-1/2)*pi and (k+1/2)*pi
def doit_bacsel(e):
   x0 = 2^(e-1)
   k0 = ceil(x0/pi+1/2)
   x1 = 2^e
   k1 = floor(x1/pi+1/2)
   t1 = RR(n((k0-1/2)*pi,200))
   if n(t1.exact_rational(),200) < n((k0-1/2)*pi,200):
      t1 = t1.nextabove()
   t1 = ZZ(t1.exact_rational()*2^(53-e))
   print ("./doit.sh 4503599627370496 " + str(t1) + " 53 " + str(e) + " 64 20")
   for k in range(k0+1,k1+1):
      t0 = t1
      t1 = RR(n((k-1/2)*pi,200))
      if n(t1.exact_rational(),200) < n((k-1/2)*pi,200):
         t1 = t1.nextabove()
      t1 = ZZ(t1.exact_rational()*2^(53-e))
      print ("./doit.sh " + str(t0) + " " + str(t1) + " 53 " + str(e) + " 64 20")
   t0 = t1
   print ("./doit.sh " + str(t0) + " 9007199254740992 53 " + str(e) + " 64 20")
   
   

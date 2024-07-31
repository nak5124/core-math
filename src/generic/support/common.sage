# SAGE file defining some common functions used elsewhere

# print like in C: 0x1.xxxp+e
def get_hex(x):
   if x < 0:
      return '-' + get_hex(-x)
   s = x.hex()
   if s[2] in ['0','1']:
      e = 0
   elif s[2] in ['2','3']:
      s = (x/2).hex()
      e = 1
   elif s[2] in ['4','5','6','7']:
      s = (x/4).hex()
      e = 2
   else:
      s = (x/8).hex()
      e = 3
   # add e to exponent
   s = s.split('p')
   e = ZZ(s[1])+e
   if e >= 0:
      s[1] = '+' + e.str()
   else:
      s[1] = e.str()
   return s[0] + 'p' + s[1]

def asuint64(x):
   assert x in RR
   s,m,e = x.sign_mantissa_exponent()
   if s==1:
      s=0
   else:
      s=1
   if abs(x) >= 2^-1022: # normal number
      m = m-2^52 # remove implicit bit
      e = e+1023+52
      return s*2^63+e*2^52+m
   # now deal with subnormal number
   m = abs(x).exact_rational()/2^-1074
   assert m in ZZ and m < 2^52
   return s*2^63+m

def fma(x,y,z):
   R = x.parent()
   assert y.parent() == R, "y.parent() == R"
   assert z.parent() == R, "z.parent() == R"
   return R(x.exact_rational()*y.exact_rational()+z.exact_rational())

exp2 = lambda x: 2^x
exp10 = lambda x: 10^x   
expm1 = lambda x: exp(x)-1
exp2m1 = lambda x: 2^x-1
exp10m1 = lambda x: 10^x-1
log10 = lambda x: log(x)/log(10)
log2 = lambda x: log(x)/log(2)
log1p = lambda x: log(1+x)
log2p1 = lambda x: log(1+x)/log(2)
log10p1 = lambda x: log(1+x)/log(10)
sinpi = lambda x: sin(pi*x)
cospi = lambda x: cos(pi*x)
tanpi = lambda x: tan(pi*x)
asinpi = lambda x: asin(x)/pi
acospi = lambda x: acos(x)/pi
atanpi = lambda x: atan(x)/pi

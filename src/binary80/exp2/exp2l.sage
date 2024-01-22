def print_T2():
   print ("static const long double T2[32][2] = {")
   R = RealField(64)
   maxerr = 0
   for i in range(2^5):
      x = n(2^(i/2^5), 200)
      h = R(x)
      l = R(x-h.exact_rational())
      print ("   {" + get_hex(h) + "L, " + get_hex(l) + "L},")
      maxerr = max(maxerr,abs(x-h.exact_rational()-l.exact_rational())/x)
   print ("};")
   print (log(maxerr)/log(2.))

def print_T1():
   print ("static const long double T1[32][2] = {")
   R = RealField(64)
   maxerr = 0
   for i in range(2^5):
      x = n(2^(i/2^10), 200)
      h = R(x)
      l = R(x-h.exact_rational())
      print ("   {" + get_hex(h) + "L, " + get_hex(l) + "L},")
      maxerr = max(maxerr,abs(x-h.exact_rational()-l.exact_rational())/x)
   print ("};")
   print (log(maxerr)/log(2.))

def print_T0():
   print ("static const long double T0[32][2] = {")
   R = RealField(64)
   maxerr = 0
   for i in range(2^5):
      x = n(2^(i/2^15), 200)
      h = R(x)
      l = R(x-h.exact_rational())
      print ("   {" + get_hex(h) + "L, " + get_hex(l) + "L},")
      maxerr = max(maxerr,abs(x-h.exact_rational()-l.exact_rational())/x)
   print ("};")
   print (log(maxerr)/log(2.))

# return the 'ulp' of the interval x, i.e., max(ulp(t)) for t in x
# this internal routine is used below
def RIFulp(x):
   return max(x.lower().ulp(),x.upper().ulp())

# analyze_P()
# err1= -116.000000000000
# err2= -98.0000000000000
# err3= -98.0000000000000
# err4= -142.000000000000
# err5= -126.999984741211
# err6= -126.000000000000
# rel. err= -83.7478216058109
# max l= 1.08420217248550444e-19
def analyze_P():
   err0 = 2^-83.748 # absolute error
   R = RealField(64)
   RI = RealIntervalField(64)
   x = RI(-2^-16,2^-16)
   p4 = RI(R("0x1.3b2ab70cf131bd7ep-7",16))
   p3 = RI(R("0x1.c6b08d6835c26dep-5",16))
   p2 = RI(R("0x1.ebfbdff82c58ea86p-3",16))
   p1 = RI(R("0x1.62e42fefa39ef358p-1",16))
   p0 = RI(1)
   # y = __builtin_fmal (p[4], x, p[3])
   y = p4*x+p3
   err1 = RIFulp(y)*x.abs().upper()^3
   print ("err1=", log(err1)/log(2.))
   # y = __builtin_fmal (y, x, p[2])
   y = y*x+p2
   err2 = RIFulp(y)*x.abs().upper()^2
   print ("err2=", log(err2)/log(2.))
   # fast_two_sum (h, l, p[1], y * x)
   z = y*x
   err3a = RIFulp(z)
   h = p1+z
   u = RIFulp(h)
   l = RI(-u,u)
   # the fast_two_sum error is bounded by 2^(1-2*p)*|h|
   err3b = 2^-127*h.abs().upper()
   err3 = (err3a+err3b)*x.abs().upper()
   print ("err3=", log(err3)/log(2.))
   # a_mul (h, &t, *h, x) [exact]
   h = h*x
   u = RIFulp(h)
   t = RI(-u,u)
   # t = __builtin_fmal (*l, x, t)
   t = l*x+t
   err4 = RIFulp(t)
   print ("err4=", log(err4)/log(2.))
   # fast_two_sum (h, l, p[0], *h)
   h = p0+h
   u = RIFulp(h)
   l = RI(-u,u)
   err5 = 2^-127*h.abs().upper()
   print ("err5=", log(err5)/log(2.))
   # *l += t;
   err6 = RIFulp(l)
   print ("err6=", log(err6)/log(2.))
   # convert err0 into relative error
   err0 = err0*(h+l).abs().upper()
   err = err0+err1+err2+err3+err4+err5+err6
   # convert into relative error
   err = err/(h+l).abs().lower()
   print ("rel. err=", log(err)/log(2.))
   print ("max l=", l.abs().upper())
   

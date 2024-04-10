# T1()
# maxerr= 3.3447210326662894e-30
# (0.982666015625000, 1.01684570312500)
def T1(xmin=0.7109375,xmax=1.421874,p=6):
   Rp = RealField(p)
   Zmin = 1
   Zmax = 1
   imin = floor((xmin-1)*128)
   imax = floor((xmax-1)*64)
   print ("imin=", imin, "imax=", imax)
   Err = 0
   print ("static const double T1[64][3] = {")
   for i in [imin..imax]:
      if i >= 0:
         xmin = RR(1 + i/64)
         xmax = RR(1 + (i+1)/64)
      else:
         xmin = RR(1 + i/128)
         xmax = RR(1 + (i+1)/128)
      rmin = Rp(1/xmax.exact_rational())
      rmax = Rp(1/xmin.exact_rational())
      r = rmin
      zbest = infinity
      while r<=rmax:
         zmin = abs(RR(r)*xmin-1)
         zmax = abs(RR(r)*xmax-1)
         z = max(zmin,zmax) # maximal distance to 1
         if (z<zbest) or (z==zbest and r==1):
            zbest = z
            rbest = r
         r = r.nextabove()
      Zmin = min(Zmin,RR(rbest)*xmin)
      Zmax = max(Zmax,RR(rbest)*xmax)
      X = rbest.exact_rational()
      # we want h multiple of 2^-38 so that e+h is exact, where
      # -16445 <= e < 16384
      e = -log(X)/log(2)
      h = RR(round(e*2^38)/2^38)
      H = h.exact_rational()
      l = RR(n(e-H,200))
      L = l.exact_rational()
      err = abs(n(H+L-e,200))
      Err = max(Err,err)
      print ("   {" + get_hex(rbest) + ", " + get_hex(h) + ", " + get_hex(l) + "}, /* i=" + str(i) + " */")
   print ("};")
   print ("maxerr=", float(err))
   return Zmin, Zmax

# T2()
# maxerr= 4.057340372917331e-29
# (0.999755859375000, 1.00024402141571)
def T2(xmin=0.982666015625000,xmax=1.016845703125,p=21):
   R = RealField(p)
   Zmin = 1
   Zmax = 1
   imin = floor((xmin-1)*2^12)
   imax = floor((xmax-1)*2^11)
   print ("imin=", imin, "imax=", imax)
   Err = 0
   print ("static const double T2[" + str(imax-imin+1) + "][3] = {")
   for i in [imin..imax]:
      if i >= 0:
         xmin = RR(1 + i/2^11)
         xmax = RR(1 + (i+1)/2^11)
      else:
         xmin = RR(1 + i/2^12)
         xmax = RR(1 + (i+1)/2^12)
      # we want |r*xmin-1| = |r*xmax-1|
      # thus 1 - r*xmin = r*xmax - 1
      # thus r = 2/(xmax + xmin)
      r = R(2 / (xmax + xmin))
      Zmin = min(Zmin,RR(r)*xmin)
      Zmax = max(Zmax,RR(r)*xmax)
      X = r.exact_rational()
      # we want h multiple of 2^-38 so that e+h is exact, where
      # -16445 <= e < 16384
      e = -log(X)/log(2)
      h = RR(round(e*2^38)/2^38)
      H = h.exact_rational()
      l = RR(n(e-H,200))
      L = l.exact_rational()
      err = abs(n(H+L-e,200))
      Err = max(Err,err)
      print ("   {" + get_hex(r) + ", " + get_hex(h) + ", " + get_hex(l) + "}, /* i=" + str(i) + " */")
   print ("};")
   print ("maxerr=", float(err))
   return Zmin, Zmax

# T1T2(6,21)
# (0.999755859375000, 1.00024402141571, 53)
def T1T2(p1,p2):
   # xh is multiple of 2^-38
   xmin, xmax = T1(p=p1)
   # r1 is multiple of 2^-p1 thus r1*xh is multiple of 2^(-38-p1)
   xmin, xmax = T2(xmin=xmin, xmax=xmax, p=p2)
   # r2 is multiple of 2^-p2 thus r2*xh is multiple of 2^(-38-p1-p2)
   e = max(1-xmin,xmax-1)
   return xmin, xmax, ceil(log(e/2^(-38-p1-p2))/log(2))
   
# return the 'ulp' of the interval x, i.e., max(ulp(t)) for t in x
# this internal routine is used below
def RIFulp(x):
   return max(x.lower().ulp(),x.upper().ulp())

def a_mul_double(a,b):
   h = a*b
   u = RIFulp(h)
   l = RIF(-u,u)
   return h, l

# analyzeP()
# err1= -129.722454159758
# err2= -113.999693586040
# err3= -113.999968548661
# err4= -89.9994919427438
# err5= -89.9999998014668
# err6= -89.0559319350096
# err7= -126.678071905112
# err8= -116.471057498903
# err9= -114.999999999999
# err10= -113.651741433865
# err= -84.3929034798149
# lmax= 1.91630753487530e-19
def analyzeP():
   err0 = 2^-84.514 # absolute error of the minimax polynomial
   p1h = RR("0x1.71547652b82fep+0",16)
   p1l = RR("0x1.777d10fa9419cp-56",16)
   p2 = RR("-0x1.71547652b82fep-1",16)
   p3 = RR("0x1.ec709dc3a03fbp-2",16)
   p4 = RR("-0x1.71547661d011fp-2",16)
   p5 = RR("0x1.2776c56eb6ca2p-2",16)
   p6 = RR("-0x1.bd761baf2853cp-3",16)
   p7 = RR("0x1.24996255a29ecp-3",16)
   p8 = RR("-0x1.6c9c75469b616p17",16)
   xh = RIF(-0.000244140625,0.00024402141571)
   xl = RIF(-2^-64,2^-64)
   x = xh+xl
   # xx = xh * xh
   xx = xh^2
   err_xx = 2*xh*xl + xl^2
   # c7 = __builtin_fma (p[8], xh, p[7])
   c7 = p8*xh+p7
   err1 = (RIFulp(c7)+(p8*xl).abs().upper())*x.abs().upper()^7
   print ("err1=", log(err1)/log(2.))
   # c5 = __builtin_fma (p[6], xh, p[5])
   c5 = p6*xh+p5
   err2 = (RIFulp(c5)+(p6*xl).abs().upper())*x.abs().upper()^5
   print ("err2=", log(err2)/log(2.))
   # c5 = __builtin_fma (c7, xx, c5)
   c5 = c7*xx+c5
   err3 = (RIFulp(c5)+(c7*err_xx).abs().upper())*x.abs().upper()^5
   print ("err3=", log(err3)/log(2.))
   # c3 = __builtin_fma (p[4], xh, p[3])
   c3 = p4*xh+p3
   err4 = (RIFulp(c3)+(p4*xl).abs().upper())*x.abs().upper()^3
   print ("err4=", log(err4)/log(2.))
   # c3 = __builtin_fma (c5, xx, c3)
   c3 = c5*xx+c3
   err5 = (RIFulp(c3)+(c5*err_xx).abs().upper())*x.abs().upper()^3
   print ("err5=", log(err5)/log(2.))
   # fast_two_sum_double (h, l, p[2], c3*xh)
   h = p2+c3*xh
   u = RIFulp(h)
   l = RIF(-u,u)
   err6a = (c3*xl).abs().upper() # neglecting xl
   err6b = h.abs().upper()*2^-105 # fast_two_sum error
   err6 = (err6a+err6b)*x.abs().upper()^2
   print ("err6=", log(err6)/log(2.))
   # d_mul_double (h, l, *h, *l, xh, xl)
   #      a_mul_double (h, &s, h_in, xh) [exact]
   #      t = fma (l_in, xh, s)
   #      l = fma (h_in, xl, t)
   h_in = h
   l_in = l
   h, s = a_mul_double (h_in, xh)
   t = l_in*xh+s
   l = h_in*xl+t
   err7 = ((l_in*xl).abs().upper()+RIFulp(t)+RIFulp(l))*x.abs().upper()
   print ("err7=", log(err7)/log(2.))
   # fast_two_sum_double (h, &t, p[0], *h)
   h = p1h+h
   u = RIFulp(h)
   t = RIF(-u,u)
   err8 = h.abs().upper()*2^-105*x.abs().upper()
   print ("err8=", log(err8)/log(2.))
   # *l += t + p[1]
   t += p1l
   l += t
   err9 = (RIFulp(t)+RIFulp(l))*x.abs().upper()
   print ("err9=", log(err9)/log(2.))
   # d_mul_double (h, l, *h, *l, xh, xl)
   #      a_mul_double (h, &s, h_in, xh) [exact]
   #      t = fma (l_in, xh, s)
   #      l = fma (h_in, xl, t)
   h_in = h
   l_in = l
   h, s = a_mul_double (h_in, xh)
   t = l_in*xh+s
   l = h_in*xl+t
   err10 = (l_in*xl).abs().upper()+RIFulp(t)+RIFulp(l)
   print ("err10=", log(err10)/log(2.))
   err = err0+err1+err2+err3+err4+err5+err6+err7+err8+err9+err10
   print ("err=", log(err)/log(2.))
   print ("lmax=", l.abs().upper())

   



         
      

#load("../../generic/support/common.sage")
# generate exact cases with y < 0
def check_exact(file=None):
   emin = -149
   emax = 127
   nsols = 0
   if file!=None:
      file = open(file,"w")
   for k in [1..10]:
      tmin = ceil(emin/2^k)
      tmax = floor(emax/2^k)
      for t in [tmin..tmax]:
         if t==0:
            continue # case x = 1
         e = t*2^k
         x = RR(2^e)
         y0 = RR(-2^-k)
         # x^y0 = 2^(-e*2^-k) = 2^(-t)
         # if y=n*y0, x^y = 2^(-n*t)
         if t<0:
            nmin = ceil(emin/-t)
            nmax = floor(emax/-t)
         else:
            nmin = ceil(emax/-t)
            nmax = floor(emin/-t)
         # we want y < 0
         nmin = max(1,nmin)
         for n in [nmin..nmax]:
            y = n*y0
            if file==None:
               print (get_hex(x), get_hex(y))
            else:
               file.write(get_hex(x) + "," + get_hex(y) + "\n")
            nsols += 1
   if file!=None:
      file.close()
   return nsols

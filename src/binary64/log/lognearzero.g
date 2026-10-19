# cr_log_fast for 1-2^-512 <= x < 1 + 2^-512 (i=511 or i=512)

@rnd = float<ieee_64,RND>;

p1 = -0x1.ffffffffffffap-2;
p2 = 0x1.555555554f4d8p-2;
p3 = -0x1.0000000537df6p-2;
p4 = 0x1.999a14758b084p-3;
p5 = -0x1.55362255e0f63p-3;

# we always have r=1, l1=l2=0

y = rnd(y_);
z = y - 1;
z2 rnd= z * z;
p45 rnd= fma(p5, z, p4);
p23 rnd= fma(p3, z, p2);
ph_1 rnd= fma(p45, z2, p23);
ph_2 rnd= fma(ph_1, z, p1);
ph_3 rnd = ph_2 * z2;

# since ee=0, we have h=z and l=ph_3

res = z+ph_3;

Z2 = z * z;
Res = z + (((p5 * z + p4) * Z2 + (p3 * z + p2)) * z + p1) * Z2;

{ y in [XMIN,XMAX] -> |(res-Res)/z| <= 3b-62 }


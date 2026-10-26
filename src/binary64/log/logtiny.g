@rnd = float<ieee_64,RND>;

c0 = -0x1p-1;
c1 = 0x1.5555555555555p-2;
c2 = -0x1.fffffffff7999p-3;
c3 = 0x1.9999999a0bbcfp-3;
c4 = -0x1.5555b5558ae3p-3;
c5 = 0x1.24915abb06a9p-3;

R2 = r * r;
r2 rnd= r*r;
P6 = c4+r*c5;
p6 rnd= c4+r*c5;
P4 = c2+r*c3;
p4 rnd= c2+r*c3;
P2 = c0+r*c1;
p2 rnd= c0+r*c1;
P4a = P4+R2*P6;
p4a rnd= p4+r2*p6;
P2a = P2+R2*P4a;
p2a rnd= p2+r2*p4a;
# fast_two_sum (&h, &l, r, r2 * p2)
T = R2*P2a;
t rnd= r2 * p2a;
h rnd= r + t;
l = rnd(-(h-(r+t)));
Y = r+T;
y = h + l;
eps = 0x1.d1p-61;
l_l rnd= l - eps*r;
l_h rnd= l + eps*r;
err_l = l_l - (l - eps*r);
err_h = l_h - (l + eps*r);

# { r in [1b-52,1b-9] -> |r2-R2| in ? } # 2^(-71.000000)
# { r in [1b-52,1b-9] -> |p6-P6| in ? } # 2^(-54.997131)
# { r in [1b-52,1b-9] -> |p4-P4| in ? } # 2^(-54.997131)
# { r in [1b-52,1b-9] -> |p2-P2| in ? } # 2^(-53.997131)
# { r in [1b-52,1b-9] -> |p4a-P4a| in ? } # 2^(-53.414062)
# { r in [1b-52,1b-9] -> |p2a-P2a| in ? } # 2^(-52.414062)
# { r in [1b-52,1b-9] \/ r in [-1b-9,-1b-53] -> t/r in [-1,1] } # fasttwosum condition
# { r in [1b-51,1b-9] -> |t/r| >= 1b-53 }
# { r in [1b-52,1b-9] -> |y-(r+t)| in ? } # 2^(-115.000000)
# { r in [1b-52,1b-9] -> |y-Y| in ? } # 2^(-69.677490)
# { r in [-1b-9,-1b-52] /\ @FIX(r,-53) -> |t/r| >= 1b-53 }
# { r in [1b-52,1b-9] \/ r in [-1b-9,-1b-53] /\ @FIX(r,-53) -> |(y-Y)/r| <= 29b-65 }
# { r in [1b-52,1b-9] \/ r in [-1b-9,-1b-53] -> |err_l/r| <= 1b-103 /\ |err_h/r| <= 1b-103 }

y-(r+t) -> l - (-(h-(r+t)));



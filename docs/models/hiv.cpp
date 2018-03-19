
$SET end = 360

$PARAM
s = 10
r = 0.03
Tmax = 1500
muT = 0.02
mub = 0.24
muV = 2.4
k1 = 2.4E-5
k2 = 3E-3
N = 1000
theta  = 1
T0 = 1000
V0 = 1E-3
NCRIT = 774


$MAIN

TAR_0 = T0;
V_0 = V0;

$CMT TAR L I V AUC

$ODE
dxdt_TAR = s - muT*TAR + r*TAR*(1-(TAR+L+I)/Tmax) - k1*V*TAR;
dxdt_L = k1*V*TAR - muT*L - k2*L;
dxdt_I = k2*L - mub*I;
dxdt_V = N*mub*I - k1*V*TAR - muV*V;
dxdt_AUC = V;
  
$TABLE
capture logV = log10(V);
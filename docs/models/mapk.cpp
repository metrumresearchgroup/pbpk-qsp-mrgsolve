$GLOBAL

//-- #include "global.h"
//--  Y = (max(0,x).^k)./(tau^k + max(0,x).^k);

double HillEQ(double x, double k, double tau) {
  double a = pow(std::max(x,0.0),k);
  return a/(pow(tau,k) + a);
}

// Created: Wed Jun 28 11:21:18 2017

$PARAM
// Created: Wed Jun 28 11:21:18 2017
// Parameters (103)
RASb   = 0    //   param 0
RASt   = 1    //   param 1
k1     = 1    //   param 2
tau1   = 1    //   param 3
k3     = 1    //   param 4
tau3   = 1    //   param 5
MEKt   = 1    //   param 6
MEKb   = 0    //   param 7
k2     = 1    //   param 8
tau2   = 1    //   param 9
ERKb   = 0    //   param 10
ERKt   = 1    //   param 11
k4     = 1    //   param 12
kFB1   = 1    //   param 13
tauFB1 = 1    //   param 14
tau4   = 1    //   param 15
S6b    = 0    //   param 16
S6t    = 1    //   param 17
k6     = 1    //   param 18
tau6   = 1    //   param 19
k5     = 1    //   param 20
tau5   = 1    //   param 21
Vmax   = 10   //   param 22
dmax   = 0    //   param 23
kg     = 1    //   param 24
taug   = 1    //   param 25
umax   = 1    //   param 26
wOR    = 0.8  //   param 27
AKTb   = 1    //   param 28
AKTt   = 1    //   param 29
k8     = 1    //   param 30
tau8   = 1    //   param 31
PI3Kb  = 1    //   param 32
PI3Kt  = 1    //   param 33
k7     = 1    //   param 34
tau7   = 1    //   param 35
ki1    = 1    //   param 36
taui1  = 1    //   param 37
ki2    = 1    //   param 38
taui2  = 1    //   param 39
kFB2   = 1    //   param 40
tauFB2 = 1    //   param 41
kFB3   = 1    //   param 42
ki4    = 1    //   param 43
tauFB3 = 1    //   param 44
taui4  = 1    //   param 45
BRAFb  = 1    //   param 46
BRAFt  = 1    //   param 47
CRAFb  = 1    //   param 48
CRAFt  = 1    //   param 49
ki3    = 1    //   param 50
taui3  = 1    //   param 51
G13    = 1    //   param 52
G23    = 1    //   param 53
RTK2b  = 1    //   param 54
RTK2t  = 1    //   param 55
ki5    = 1    //   param 56
taui5  = 1    //   param 57
G14    = 1    //   param 58
RTK1b  = 1    //   param 59
RTK1t  = 1    //   param 60
kFB4   = 1    //   param 61
tauFB4 = 1    //   param 62
G24    = 1    //   param 63
F1     = 1    //   param 64
V1     = 1    //   param 65
ka1    = 1    //   param 66
ke1    = 1    //   param 67
p1     = 1    //   param 68
G33    = 1    //   param 69
G34    = 1    //   param 70
RTK3b  = 1    //   param 71
F2     = 1    //   param 72
V2     = 1    //   param 73
ka2    = 1    //   param 74
ke2    = 1    //   param 75
F3     = 1    //   param 76
V3     = 1    //   param 77
ka3    = 1    //   param 78
ke3    = 1    //   param 79
p2     = 1    //   param 80
p3     = 1    //   param 81
F4     = 1    //   param 82
V4     = 1    //   param 83
ka4    = 1    //   param 84
ke4    = 1    //   param 85
p4     = 1    //   param 86
wRAS   = 1    //   param 87
RTK3t  = 1    //   param 88
F5     = 1    //   param 89
V5     = 1    //   param 90
ka5    = 1    //   param 91
ke5    = 1    //   param 92
p5     = 1    //   param 93
r4     = 1    //   param 94
r2     = 1    //   param 95
r3     = 1    //   param 96
r1     = 1    //   param 97
q2     = 1    //   param 98
V3b    = 1    //   param 99
r5     = 1    //   param 100
Gdusp  = 1    //   param 101
Gspry  = 1    //   param 102

$INIT
// Created: Wed Jun 28 11:21:18 2017
// Initial conditions (17)
TD1         = 0  //   species 6
CELLS       = 1  //   species 7
FB1         = 0  //   species 8
FB2         = 0  //   species 9
FB3         = 0  //   species 10
FB4         = 0  //   species 20
RTK1i_blood = 0  //   species 22
RAFi_gut    = 0  //   species 23
RAFi_blood  = 0  //   species 24
MEKi_gut    = 0  //   species 25
MEKi_blood  = 0  //   species 26
ERKi_gut    = 0  //   species 27
ERKi_blood  = 0  //   species 28
AKTi_gut    = 0  //   species 29
AKTi_blood  = 0  //   species 30
MEKi_V3     = 0  //   species 31
RTK1i_gut   = 0  //   species 32

$ODE

// Created: Wed Jun 28 11:21:18 2017
// RULES (21)
double RTK1i_C = RTK1i_blood / V1;
double RAFi_C  = RAFi_blood / V2;
double MEKi_C  = MEKi_blood / V3;
double ERKi_C  = ERKi_blood / V4;
double AKTi_C  = AKTi_blood / V5;
double RTK1i   = p1 * RTK1i_C;
double RAFi    = p2 * RAFi_C;
double MEKi    = p3 * MEKi_C;
double ERKi    = p4 * ERKi_C;
double AKTi    = p5 * AKTi_C;
double RTK1    = (RTK1b + (RTK1t - RTK1b) * (1 - G13 * HillEQ(FB3, kFB3, tauFB3)) * (1 - G14 * HillEQ(FB4, kFB4, tauFB4))) * (1 - HillEQ(RTK1i, ki1, taui1));
double RTK2    = RTK2b + (RTK2t - RTK2b) * (1 - G23 * HillEQ(FB3, kFB3, tauFB3)) * (1 - G24 * HillEQ(FB4, kFB4, tauFB4));
double RTK3    = RTK3b + (RTK3t - RTK3b) * (1 - G33 * HillEQ(FB3, kFB3, tauFB3)) * (1 - G34 * HillEQ(FB4, kFB4, tauFB4));
double RAS     = (RASb + (RASt - RASb) * HillEQ(RTK1 + RTK2, k1, tau1)) * (1 - Gspry * HillEQ(FB2, kFB2, tauFB2));
double BRAF    = (BRAFb + (BRAFt - BRAFb) * HillEQ(RAS, k2, tau2)) * (1 - HillEQ(RAFi, ki2, taui2));
double CRAF    = CRAFb + (CRAFt - CRAFb) * HillEQ(RAS, k5, tau5);
double MEK     = (MEKb + (MEKt - MEKb) * HillEQ(BRAF + CRAF, k3, tau3)) * (1 - HillEQ(MEKi, ki3, taui3));
double ERK     = (ERKb + (ERKt - ERKb) * HillEQ(MEK, k4, tau4)) * (1 - Gdusp * HillEQ(FB1, kFB1, tauFB1)) * (1 - HillEQ(ERKi, ki4, taui4));
double PI3K    = PI3Kb + (PI3Kt - PI3Kb) * HillEQ(RTK3 + wRAS * RAS, k7, tau7);
double AKT     = (AKTb + (AKTt - AKTb) * HillEQ(PI3K, k8, tau8)) * (1 - HillEQ(AKTi, ki5, taui5));
double S6      = S6b + (S6t - S6b) * HillEQ(wOR * ERK + (1 - wOR) * AKT, k6, tau6);

// Created: Wed Jun 28 11:21:18 2017
// Reactions (22)
double fb1        = r1 * (ERK - FB1);
double fb2        = r2 * (ERK - FB2);
double fb3        = r3 * (ERK - FB3);
double TD         = r5 * (S6 - TD1);
double cell       = (umax * HillEQ(TD1, kg, taug) - dmax) * CELLS * (1 - CELLS / Vmax);
double fb4        = r4 * (AKT - FB4);
double PK1a_RAFi  = ka2 * F2 * RAFi_gut;
double PK1a_MEKi  = ka3 * F3 * MEKi_gut;
double PK1a_ERKi  = ka4 * F4 * ERKi_gut;
double PK1a_AKTi  = ka5 * F5 * AKTi_gut;
double PK2_RTK1i  = ke1 * RTK1i_blood;
double PK2_RAFi   = ke2 * RAFi_blood;
double PK2_MEKi   = ke3 * MEKi_blood;
double PK2_ERKi   = ke4 * ERKi_blood;
double PK2_AKTi   = ke5 * AKTi_blood;
double PK3_MEKi   = q2 / V3 * MEKi_blood - q2 / V3b * MEKi_V3;
double PK1a_RTK1i = ka1 * F1 * RTK1i_gut;
double PK1b_RTK1i = ka1 * (1 - F1) * RTK1i_gut;
double PK1b_RAFi  = ka2 * (1 - F2) * RAFi_gut;
double PK1b_MEKi  = ka3 * (1 - F3) * MEKi_gut;
double PK1b_ERKi  = ka4 * (1 - F4) * ERKi_gut;
double PK1b_AKTi  = ka5 * (1 - F5) * AKTi_gut;


// Created: Wed Jun 28 11:21:18 2017
// ODEs (17)
dxdt_AKTi_blood  =   PK1a_AKTi - PK2_AKTi;
dxdt_AKTi_gut    = - PK1a_AKTi - PK1b_AKTi;
dxdt_CELLS       =   cell;
dxdt_ERKi_blood  =   PK1a_ERKi - PK2_ERKi;
dxdt_ERKi_gut    = - PK1a_ERKi - PK1b_ERKi;
dxdt_FB1         =   fb1;
dxdt_FB2         =   fb2;
dxdt_FB3         =   fb3;
dxdt_FB4         =   fb4;
dxdt_MEKi_blood  =   PK1a_MEKi - PK2_MEKi - PK3_MEKi;
dxdt_MEKi_gut    = - PK1a_MEKi - PK1b_MEKi;
dxdt_MEKi_V3     =   PK3_MEKi;
dxdt_RAFi_blood  =   PK1a_RAFi - PK2_RAFi;
dxdt_RAFi_gut    = - PK1a_RAFi - PK1b_RAFi;
dxdt_RTK1i_blood =   PK1a_RTK1i - PK2_RTK1i;
dxdt_RTK1i_gut   = - PK1a_RTK1i - PK1b_RTK1i;
dxdt_TD1         =   TD;

$TABLE
capture TUMOR = CELLS;
capture GDC = ERKi;

$CAPTURE ERKi ERKi_C


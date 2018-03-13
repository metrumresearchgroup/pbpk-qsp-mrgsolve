[PROB]
# Yoshikado et al. (2016)

https://www.ncbi.nlm.nih.gov/pubmed/27170342

- Title: __Quantitative Analyses of Hepatic OATP-Mediated
Interactions Between Statins and Inhibitors Using PBPK
Modeling With a Parameter Optimiaztion Method__
- Reference: CP\&T vol. 100 no. 5 pp. 513-23 11/2016
- Parameters: 40
- Compartments: 31

[CMT] 

// Dosing
gut igut

// Statin compartments
cent mus adi ski
ehc1 ehc2 ehc3
he1 he2 he3 he4 he5
hc1 hc2 hc3 hc4 hc5

// CsA compartments
icent
me se ae
mc sc ac
iliv1 iliv2 iliv3 iliv4 iliv5

[PARAM] // CSA
iKp_mus = 2.98
iKp_adi = 17.3
iKp_ski = 13.6
iKp_liv = 16.7
ifb = 0.06
ikiu = 0.0118
imw = 1202.61

PSmus = 245/70
PSski = 37.4/70
PSadi = 10.2/70
ifhCLint = 0.587/70

ifafg = 0.572
iClr = 0
ika = 0.999
itlag =  0.254 


Kp_ski = 0.481
Kp_mus = 0.113
Kp_adi = 0.086
CLr = 0.0
Vcent = 0.075
fafg = 1.0
ktr = 0.679
ka = 1.06
fb = 0.008
fh = 0.035
fbCLintall = 51.6/70
fbile = 0.330
gamma = 0.244
beta = 0.8
Rdiff = 0.0345
tlag = 1

[PARAM]
Qh   = 1.200
Qmus = 0.642
Qski = 0.257
Qadi = 0.223

Vliv = 0.0241
Vmus = 0.4290
Vski = 0.1110
Vadi = 0.1430

exFliv = 0.278
exFmus = 0.146
exFski = 0.321
exFadi = 0.145


[MAIN]

if(NEWIND <=1) {
  double CLintall = fbCLintall/fb;
  double PSact = 1.0/(1.0+Rdiff)*CLintall/beta;
  double PSdiffi = Rdiff/(1.0+Rdiff)*CLintall/beta;
  double PSdiffe = Rdiff/(1.0+Rdiff)/gamma*CLintall/beta;
  double CLint = CLintall/(1.0-beta)*Rdiff/(1.0+Rdiff)/gamma;
  
  double Vme = Vmus*exFmus;
  double Vae = Vadi*exFadi;
  double Vse = Vski*exFski;
  double Vmc = Vmus-Vme;
  double Vac = Vadi-Vae;
  double Vsc = Vski-Vse;
  double dVliv = Vliv/5.0;
  double ikitot = imw*ikiu/ifb;
}

// ALAG_gut = tlag;
//ALAG_igut = itlag;

[ODE]

// Statin concentrations
double Ccent = cent/Vcent;
double Cmus  = mus/Vmus;
double Cski  = ski/Vski;
double Cadi  = adi/Vadi;

// Volume liv
double Vhe = dVliv*exFliv;
double Vhc = dVliv*(1-exFliv);

double Chc1 = hc1/Vhc;
double Chc2 = hc2/Vhc;
double Chc3 = hc3/Vhc;
double Chc4 = hc4/Vhc;
double Chc5 = hc5/Vhc;

double Che1 = he1/Vhe;
double Che2 = he2/Vhe;
double Che3 = he3/Vhe;
double Che4 = he4/Vhe;
double Che5 = he5/Vhe;

// CsA concentrations
double iCcent = icent/Vcent;

double Cme = me/Vme;
double Cse = se/Vse;
double Cae = ae/Vae;

double Cmc = mc/Vmc;
double Csc = sc/Vsc;
double Cac = ac/Vac;

double iCliv1 = iliv1/dVliv;
double iCliv2 = iliv2/dVliv;
double iCliv3 = iliv3/dVliv;
double iCliv4 = iliv4/dVliv;
double iCliv5 = iliv5/dVliv;

dxdt_igut = -ika/ifafg*igut;

dxdt_icent = 
  Qh*iCliv5/iKp_liv 
  - Qh*iCcent 
  - iClr*iCcent 
  - Qmus*(iCcent-Cme) 
  - Qski*(iCcent-Cse) 
  - Qadi*(iCcent-Cae);
  
dxdt_me = Qmus*(iCcent-Cme) - PSmus*ifb*(Cme-Cmc/iKp_mus);
dxdt_se = Qski*(iCcent-Cse) - PSski*ifb*(Cse-Csc/iKp_ski);
dxdt_ae = Qadi*(iCcent-Cae) - PSadi*ifb*(Cae-Cac/iKp_adi);
  
dxdt_mc = PSmus*ifb*(Cme-Cmc/iKp_mus);
dxdt_sc = PSski*ifb*(Cse-Csc/iKp_ski);
dxdt_ac = PSadi*ifb*(Cae-Cac/iKp_adi);
  
dxdt_iliv1 = Qh*(iCcent-iCliv1/iKp_liv) - (ifhCLint/5.0)*iCliv1 + ika*igut;
dxdt_iliv2 = Qh*(iCliv1-iCliv2)/iKp_liv - (ifhCLint/5.0)*iCliv2;
dxdt_iliv3 = Qh*(iCliv2-iCliv3)/iKp_liv - (ifhCLint/5.0)*iCliv3;
dxdt_iliv4 = Qh*(iCliv3-iCliv4)/iKp_liv - (ifhCLint/5.0)*iCliv4;
dxdt_iliv5 = Qh*(iCliv4-iCliv5)/iKp_liv - (ifhCLint/5.0)*iCliv5;
  
// CsA effect on Statin
double csai1 = 1.0+(iCliv1/iKp_liv)/ikitot;
double csai2 = 1.0+(iCliv2/iKp_liv)/ikitot;
double csai3 = 1.0+(iCliv3/iKp_liv)/ikitot;
double csai4 = 1.0+(iCliv4/iKp_liv)/ikitot;
double csai5 = 1.0+(iCliv5/iKp_liv)/ikitot;


double hex2 = fh*(PSdiffe/5.0);
dxdt_he1 = Qh*(Ccent-Che1)-(fb*(PSact/csai1+PSdiffi)/5.0)*Che1+hex2*Chc1 + ka*gut;
dxdt_he2 = Qh*(Che1 -Che2)-(fb*(PSact/csai2+PSdiffi)/5.0)*Che2+hex2*Chc2;
dxdt_he3 = Qh*(Che2 -Che3)-(fb*(PSact/csai3+PSdiffi)/5.0)*Che3+hex2*Chc3;
dxdt_he4 = Qh*(Che3 -Che4)-(fb*(PSact/csai4+PSdiffi)/5.0)*Che4+hex2*Chc4;
dxdt_he5 = Qh*(Che4 -Che5)-(fb*(PSact/csai5+PSdiffi)/5.0)*Che5+hex2*Chc5;


double hcx2 = fh*((PSdiffe+CLint)/5.0);
dxdt_hc1 = fb*((PSact/csai1+PSdiffi)/5.0)*Che1 - hcx2*Chc1;
dxdt_hc2 = fb*((PSact/csai2+PSdiffi)/5.0)*Che2 - hcx2*Chc2;
dxdt_hc3 = fb*((PSact/csai3+PSdiffi)/5.0)*Che3 - hcx2*Chc3;
dxdt_hc4 = fb*((PSact/csai4+PSdiffi)/5.0)*Che4 - hcx2*Chc4;
dxdt_hc5 = fb*((PSact/csai5+PSdiffi)/5.0)*Che5 - hcx2*Chc5;

dxdt_cent = 
  Qh*Che5 
  - Qh*Ccent 
  - CLr*Ccent 
  - Qmus*(Ccent-Cmus/Kp_mus) 
  - Qski*(Ccent-Cski/Kp_ski) 
  - Qadi*(Ccent-Cadi/Kp_adi);

dxdt_mus = Qmus*(Ccent-Cmus/Kp_mus);
dxdt_ski = Qski*(Ccent-Cski/Kp_ski);
dxdt_adi = Qadi*(Ccent-Cadi/Kp_adi);


dxdt_gut  = ktr*ehc3 - ka/fafg*gut;
dxdt_ehc1 = fbile*fh*(CLint/5.0)*(Chc1+Chc2+Chc3+Chc4+Chc5)-ktr*ehc1;
dxdt_ehc2 = ktr*(ehc1-ehc2);
dxdt_ehc3 = ktr*(ehc2-ehc3);



[TABLE]
capture CP =  cent/Vcent;
capture CSA = icent/Vcent;
capture CSAliv = iliv1/dVliv;


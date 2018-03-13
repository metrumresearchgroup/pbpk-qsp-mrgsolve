[ PROB ]

1: Asaumi R, Toshimoto K, Tobe Y, Hashizume K, Nunoya KI, Imawaka H, Lee W,
Sugiyama Y. Comprehensive PBPK Model of Rifampicin for Quantitative Prediction of
Complex Drug-Drug Interactions: CYP3A/2C9 Induction and OATP Inhibition Effects. 
CPT Pharmacometrics Syst Pharmacol. 2018 Jan 25. doi: 10.1002/psp4.12275. [Epub
ahead of print] PubMed PMID: 29368402.

https://www.ncbi.nlm.nih.gov/pubmed/29368402
http://onlinelibrary.wiley.com/doi/10.1002/psp4.12275/abstract

[ PARAM ]

Rdif  = 0.129
beta  = 0.2    // 0.2/0.5/0.8
gamma = 0.778

Km_u_uptake = 0.146 // 1.23 // ug/mL

SFKp  = 6.65  // 1 
mSFKp = 0.201 // 1

Emax_UGT_RIF   = 1.34   // 1.00
EC50_u_UGT_RIF = 0.0526 // ug/mL

kdeg_UGT_liver = 0.0158 // per h
kdeg_UGT_ent   = 0.0288 // per h

fm_UGT_liver = 0.759 // same
fm_UGT_ent   = 0.759 // same

Emax_CYP3A4_RIF   = 4.5700  // 12.3
EC50_u_CYP3A4_RIF = 0.0526 // ug/mL

kdeg_CYP3A4_liver = 0.0158 // per h
kdeg_CYP3A4_ent   = 0.0288 // per hr

fm_CYP3A4_liver = 0.93
fm_CYP3A4_ent   = 1.00

mCLperm_gut_kg = 0.151


[ PARAM ]
fB  = 0.0778
mfB = 0.0545
fH  = 0.0814
fE  = 0.115

Fa  = 1.000
mFa = 1.000
Fg  = 0.943

[ PARAM ]
Kp_skin    = 0.326
Kp_muscle  = 0.0947
Kp_adipose = 0.0629
Kp_serosa  = 0.200

mKp_liver   = 6.96
mKp_muscle  = 4.00
mKp_skin    = 20.4
mKp_adipose = 34.4

[PARAM]
fBCLint_all_kg = 0.251 // 0.204 // L/h/kg
mfBCLint_kg    = 0.469 // 0.528 // L/h/kg
mfECLint_E_kg  = 0.107 // L/h/kg

PSdif_E_kg = 0.161 // 0.143 // L/h/kg
CLrenal_kg = 0.011 // L/h/kg
mCLrenal   = 0.000

[ PARAM ]
Qvilli_kg   = 0.257 // L/h/kg
Qh_kg       = 1.240 // L/h/kg
Qmuscle_kg  = 0.642 // L/h/kg
Qskin_kg    = 0.257 // L/h/kg
Qadipose_kg = 0.223 // L/h/kg
Qserosa_kg  = 0.274 // L/h/kg
Qportal_kg  = 0.531 // L/h/kg

[ PARAM ]
VHE_kg       = 0.0067  // L/kg (Vi)
VHC_kg       = 0.0174  // L/kg (Vh)
Vcentral_kg  = 0.0743  // L/kg (VbRif)
mVcentral_kg = 0.571   // L/kg (estimated)
Vskin_kg     = 0.111   // L/kg
Vadipose_kg  = 0.143   // L/kg
Vmuscle_kg   = 0.429   // L/kg
Vserosa_kg   = 0.00893 // L/kg
Vent_kg      = 0.00739 // L/kg
Vmucblood_kg = 0.00099 // L/kg
Vportal_kg   = 0.001   // L/kg

[PARAM]
ka  = 37.6  // 3.26 // per hr
mka = 1.29  // 5.51 per hr
WT  = 80    // Not sure 

[ TABLE ] 
capture Cmidazolam = 1000*mCcentral;capture Rflux  = 

[CAPTURE] Ccentral mCcentral

[ MAIN ]
if(NEWIND <= 1) {
  // -------------------------------------------
  double fBCLint_all = fBCLint_all_kg*WT;
  double CLint_all   = fBCLint_all / fB;
  double mfBCLint    = mfBCLint_kg*WT;
  double mCLperm_gut = mCLperm_gut_kg*WT;
  double mfECLint_E  = mfECLint_E_kg*WT;
  double CLrenal     = CLrenal_kg*WT;
  double PSdif_E     = PSdif_E_kg*WT;
  // -------------------------------------------
  double Qvilli    = Qvilli_kg*WT;
  double Qh        = Qh_kg*WT;
  double Qmuscle   = Qmuscle_kg*WT;
  double Qskin     = Qskin_kg*WT;
  double Qadipose  = Qadipose_kg*WT;
  double Qserosa   = Qserosa_kg*WT;
  double Qhart     = Qh - Qserosa - Qvilli;
  double Qportal   = Qportal_kg*WT;
  // -------------------------------------------
  double VHE = VHE_kg*WT;
  double VHC = VHC_kg*WT;
  double Vcentral = Vcentral_kg*WT;
  double Vskin = Vskin_kg*WT;
  double Vadipose = Vadipose_kg*WT;
  double Vmuscle = Vmuscle_kg*WT;
  double Vserosa = Vserosa_kg*WT;
  double Vent = Vent_kg*WT;
  double Vmucblood = Vmucblood_kg*WT;
  double mVcentral = mVcentral_kg*WT;
  double Vportal = Vportal_kg*WT;
  // -------------------------------------------
  double Vmax_uptake = 1.0 / (1 + Rdif) * CLint_all / beta * Km_u_uptake;
  double PSdif_inf = Rdif /  (1 + Rdif) * CLint_all / beta;
  double PSdif_eff = Rdif /  (1 + Rdif) * CLint_all / beta / gamma;
  double CLint = Rdif / (1 + Rdif) * CLint_all / (1 - beta) / gamma;
  double Qgut = fE * PSdif_E * Qvilli / (Qvilli + fB * PSdif_E);
  double mQgut = Qvilli * mCLperm_gut / (Qvilli + mCLperm_gut);
  double CLint_E = (Qgut * (1.0 / Fg - 1.0) - (1.0- Fa) * fE * PSdif_E * 20.0) / fE;
}

[CMT]
Xgutlumen Mgutlumen
central Cmuscle Cskin Cadipose Cserosa Cmucblood Cent
CHE1 CHE2 CHE3 CHE4 CHE5
CHC1 CHC2 CHC3 CHC4 CHC5

mcentral mCmuscle mCskin mCadipose
CLIV1 CLIV2 CLIV3 CLIV4 CLIV5  
Cportal 
  
[INIT]
UGT_ratio_HC1 = 1
UGT_ratio_HC2 = 1
UGT_ratio_HC3 = 1
UGT_ratio_HC4 = 1
UGT_ratio_HC5 = 1
UGT_ratio_ent = 1

CYP3A4_ratio_HC1 = 1
CYP3A4_ratio_HC2 = 1
CYP3A4_ratio_HC3 = 1
CYP3A4_ratio_HC4 = 1
CYP3A4_ratio_HC5 = 1
CYP3A4_ratio_ent = 1

[ ODE ]

double Ccentral = central/Vcentral;

dxdt_central = 
  Qh       * CHE5 - 
  Qhart    * Ccentral - 
  Qserosa  * Ccentral - 
  Qvilli   * Ccentral - 
  CLrenal  * Ccentral + 
  Qmuscle  * (Cmuscle  / (SFKp * Kp_muscle)  - Ccentral) + 
  Qskin    * (Cskin    / (SFKp * Kp_skin)    - Ccentral) + 
  Qadipose * (Cadipose / (SFKp * Kp_adipose) - Ccentral);

dxdt_Cmuscle = 
  (1.0/Vmuscle) * Qmuscle * (Ccentral - Cmuscle / (SFKp * Kp_muscle));

dxdt_Cskin = 
  (1.0/Vskin) * Qskin * (Ccentral - Cskin / (SFKp * Kp_skin));

dxdt_Cadipose = 
  (1.0/Vadipose) * Qadipose * (Ccentral - Cadipose / (SFKp * Kp_adipose));

dxdt_Cserosa = 
  (1.0/Vserosa) * Qserosa * (Ccentral - Cserosa / (SFKp * Kp_serosa));

dxdt_Cmucblood = 
  Qvilli * (Ccentral - Cmucblood) + 
  fE * PSdif_E * Cent - fB * PSdif_E * Cmucblood;

dxdt_Cmucblood = dxdt_Cmucblood * (1/Vmucblood);

dxdt_Xgutlumen = 
  - ka / Fa * Xgutlumen + fE * PSdif_E * 20 * Cent;
  
dxdt_Cent = 
  ka * Xgutlumen + 
  fB * PSdif_E * Cmucblood - 
  fE * (PSdif_E * 21 + 
  CLint_E * (1 + fm_UGT_ent * (UGT_ratio_ent - 1))) * Cent;

dxdt_Cent = dxdt_Cent * (1/Vent);

dxdt_UGT_ratio_HC1 = 
  kdeg_UGT_liver * 
  (1 + Emax_UGT_RIF * fH * CHC1 / (fH * CHC1 + EC50_u_UGT_RIF) - UGT_ratio_HC1);

dxdt_UGT_ratio_HC2 = 
  kdeg_UGT_liver * 
  (1 + Emax_UGT_RIF * fH * CHC2 / (fH * CHC2 + EC50_u_UGT_RIF) - UGT_ratio_HC2);

dxdt_UGT_ratio_HC3 = 
  kdeg_UGT_liver *
  (1 + Emax_UGT_RIF * fH * CHC3 / (fH * CHC3 + EC50_u_UGT_RIF) - UGT_ratio_HC3);

dxdt_UGT_ratio_HC4 = 
  kdeg_UGT_liver * 
  (1 + Emax_UGT_RIF * fH * CHC4 / (fH * CHC4 + EC50_u_UGT_RIF) - UGT_ratio_HC4);

dxdt_UGT_ratio_HC5 = 
  kdeg_UGT_liver * 
  (1 + Emax_UGT_RIF * fH * CHC5 / (fH * CHC5 + EC50_u_UGT_RIF) - UGT_ratio_HC5);

dxdt_UGT_ratio_ent = 
  kdeg_UGT_ent * 
  (1 + Emax_UGT_RIF * fE * Cent / (fE * Cent + EC50_u_UGT_RIF) - UGT_ratio_ent);

dxdt_CHE1 = 
  Qhart  * Ccentral + 
  Qvilli * Cmucblood + 
  Qserosa * Cserosa / (SFKp * Kp_serosa) - 
  Qh * CHE1 + 
  (fH * PSdif_eff * CHC1 - 
   fB * (Vmax_uptake / (Km_u_uptake + fB * CHE1) + PSdif_inf) * CHE1) / 5.0;   

dxdt_CHE1 = dxdt_CHE1 * (5.0/VHE);

dxdt_CHE2 = 
  Qh * (CHE1 - CHE2) + 
  (fH * PSdif_eff * CHC2 - 
   fB * (Vmax_uptake / (Km_u_uptake + fB * CHE2) + PSdif_inf) * CHE2) / 5.0;

dxdt_CHE2 = dxdt_CHE2 * (5.0/VHE);

dxdt_CHE3 = 
  Qh * (CHE2 - CHE3) + 
  (fH * PSdif_eff * CHC3 - 
   fB * (Vmax_uptake / (Km_u_uptake + fB * CHE3) + PSdif_inf) * CHE3) / 5.0;

dxdt_CHE3 = dxdt_CHE3 * (5.0/VHE);

dxdt_CHE4 = 
  Qh * (CHE3 - CHE4) + 
  (fH * PSdif_eff * CHC4 - 
   fB * (Vmax_uptake / (Km_u_uptake + fB * CHE4) + PSdif_inf) * CHE4) / 5.0;

dxdt_CHE4 = dxdt_CHE4 * (5.0/VHE);

dxdt_CHE5 = 
  Qh * (CHE4 - CHE5) + 
  (fH * PSdif_eff * CHC5 - 
   fB * (Vmax_uptake / (Km_u_uptake + fB * CHE5) + PSdif_inf) * CHE5) / 5.0;// (i = 2~5)

dxdt_CHE5 = dxdt_CHE5 * (5.0/VHE);

dxdt_CHC1 = 
  (5.0/VHC) *
  (fB * (Vmax_uptake / (Km_u_uptake + fB * CHE1) + PSdif_inf) * CHE1 - 
   fH * PSdif_eff * CHC1 - 
   fH * CLint * (1 + fm_UGT_liver * (UGT_ratio_HC1 - 1)) * CHC1) / 5.0;

dxdt_CHC2 = 
  (5.0/VHC) *
  (fB * (Vmax_uptake / (Km_u_uptake + fB * CHE2) + PSdif_inf) * CHE2 - 
   fH * PSdif_eff * CHC2 - 
   fH * CLint * (1 + fm_UGT_liver * (UGT_ratio_HC2 - 1)) * CHC2) / 5.0;

dxdt_CHC3 = 
  (5.0/VHC) *
  (fB * (Vmax_uptake / (Km_u_uptake + fB * CHE3) + PSdif_inf) * CHE3 - 
   fH * PSdif_eff * CHC3 - 
   fH * CLint * (1 + fm_UGT_liver * (UGT_ratio_HC3 - 1)) * CHC3) / 5.0;

dxdt_CHC4 = 
  (5.0/VHC) * 
  (fB * (Vmax_uptake / (Km_u_uptake + fB * CHE4) + PSdif_inf) * CHE4 - 
   fH * PSdif_eff * CHC4 -
   fH * CLint * (1 + fm_UGT_liver * (UGT_ratio_HC4 - 1)) * CHC4) / 5.0;

dxdt_CHC5 = 
  (5.0/VHC) *
  (fB * (Vmax_uptake / (Km_u_uptake + fB * CHE5) + PSdif_inf) * CHE5 - 
   fH * PSdif_eff * CHC5 - 
   fH * CLint * (1 + fm_UGT_liver * (UGT_ratio_HC5 - 1)) * CHC5) / 5.0;

double mCcentral = mcentral/mVcentral;
dxdt_mcentral = 
  Qh * (CLIV5 / (mSFKp * mKp_liver)) - 
  (Qh-Qportal) * mCcentral +
  Qmuscle      * (mCmuscle  / (mSFKp * mKp_muscle)  - mCcentral) + 
  Qskin        * (mCskin    / (mSFKp * mKp_skin)    - mCcentral) +
  Qadipose     * (mCadipose / (mSFKp * mKp_adipose) - mCcentral) -
  Qportal      * mCcentral - 
  mCLrenal     * mCcentral;

dxdt_CYP3A4_ratio_HC1 =
  kdeg_CYP3A4_liver *
  (1 + Emax_CYP3A4_RIF * fH * CHC1 / (fH * CHC1 + EC50_u_CYP3A4_RIF) - CYP3A4_ratio_HC1);
dxdt_CYP3A4_ratio_HC2 =
  kdeg_CYP3A4_liver *
  (1 + Emax_CYP3A4_RIF * fH * CHC2 / (fH * CHC2 + EC50_u_CYP3A4_RIF) - CYP3A4_ratio_HC2);
dxdt_CYP3A4_ratio_HC3 =
  kdeg_CYP3A4_liver *
  (1 + Emax_CYP3A4_RIF * fH * CHC3 / (fH * CHC3 + EC50_u_CYP3A4_RIF) - CYP3A4_ratio_HC3);
dxdt_CYP3A4_ratio_HC4 =
  kdeg_CYP3A4_liver *
  (1 + Emax_CYP3A4_RIF * fH * CHC4 / (fH * CHC4 + EC50_u_CYP3A4_RIF) - CYP3A4_ratio_HC4);
dxdt_CYP3A4_ratio_HC5 =
  kdeg_CYP3A4_liver *
  (1 + Emax_CYP3A4_RIF * fH * CHC5 / (fH * CHC5 + EC50_u_CYP3A4_RIF) - CYP3A4_ratio_HC5);
dxdt_CYP3A4_ratio_ent = 
  kdeg_CYP3A4_ent * 
  (1 + Emax_CYP3A4_RIF * fE * Cent / (fE * Cent + EC50_u_CYP3A4_RIF) - CYP3A4_ratio_ent);
          
dxdt_CLIV1 = 
  (Qh-Qportal) * mCcentral + 
  Qportal * Cportal - 
  Qh * CLIV1 / (mSFKp * mKp_liver) - 
  mfBCLint * (1 + fm_CYP3A4_liver * (CYP3A4_ratio_HC1 - 1)) / 5 * CLIV1 / 
  (mSFKp * mKp_liver);

dxdt_CLIV1 = dxdt_CLIV1 * (5/(VHE+VHC));
  
dxdt_CLIV2 = 
  (Qh * (CLIV1 - CLIV2) - 
   mfBCLint * (1 + fm_CYP3A4_liver * (CYP3A4_ratio_HC2 - 1)) / 5 * CLIV2) / 
   (mSFKp * mKp_liver); 
  
dxdt_CLIV2 = dxdt_CLIV2 * (5/(VHE+VHC));   
 
dxdt_CLIV3 = 
  (Qh * (CLIV2 - CLIV3) - 
  mfBCLint * (1 + fm_CYP3A4_liver * (CYP3A4_ratio_HC3 - 1)) / 5 * CLIV3) / 
  (mSFKp * mKp_liver); 
 
dxdt_CLIV3 = dxdt_CLIV3 * (5/(VHE+VHC));   
 
dxdt_CLIV4 = 
  (Qh * (CLIV3 - CLIV4) - 
  mfBCLint * (1 + fm_CYP3A4_liver * (CYP3A4_ratio_HC4 - 1)) / 5 * CLIV4) / 
  (mSFKp * mKp_liver); 
 
dxdt_CLIV4 = dxdt_CLIV4 * (5/(VHE+VHC));    

dxdt_CLIV5 = 
  (Qh * (CLIV4 - CLIV5) - 
  mfBCLint * (1 + fm_CYP3A4_liver * (CYP3A4_ratio_HC5 - 1)) / 5 * CLIV5) / 
  (mSFKp * mKp_liver); 

dxdt_CLIV5 = dxdt_CLIV5 * (5/(VHE+VHC));    

dxdt_Cportal = 
  Qportal * (mCcentral - Cportal) + 
  mka * mQgut / (mQgut + mfECLint_E * (1 + fm_CYP3A4_ent * (CYP3A4_ratio_ent - 1))) * Mgutlumen; 
dxdt_Cportal = dxdt_Cportal * (1/Vportal); 
  
dxdt_mCmuscle = 
  (1/Vmuscle) * Qmuscle * (mCcentral - mCmuscle / (mSFKp * mKp_muscle));

dxdt_mCskin = 
  (1/Vskin) * Qskin * (mCcentral - mCskin / (mSFKp * mKp_skin));

dxdt_mCadipose = 
  (1/Vadipose) * Qadipose * (mCcentral - mCadipose / (mSFKp * mKp_adipose));

dxdt_Mgutlumen = -mka/mFa * Mgutlumen;


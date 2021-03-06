#parameter values of case report model

#kinetic parameters
gamma_EGFR=0.02
gamma_IGF1R=0.02
gamma_ERBB4=0.02
gamma_cKIT=0.02
kf_PI3K_active=0.005
k_p90Rsk_Erk=0.0213697
KM_p90Rsk_Erk=763523
k_SOS_E=694.731
KM_SOS_E=6086070
k_SOS_ERBB4=7
KM_SOS_ERBB4=400000
k_SOS_cKIT=7
KM_SOS_cKIT=400000
k_Ras_SOS=32.344
KM_Ras_SOS=35954.3
k_ERK_MEK=9.85367
KM_ERK_MEK=1007340
k_DSOS_p90Rsk=161197
KM_DSOS_p90Rsk=896896
k_SOS_I=500
KM_SOS_I=100000
k_PIK3_IGF1R=10.6737
KM_PIK3_IGF1R=184912
k_PIK3_EGFR=10.6737
KM_PIK3_EGFR=184912
k_Akt_PIK3=0.0566279
KM_Akt_PIK3=653951
kd_Akt=0.005
k_Erk_PP2A=8.8912
KM_Erk_PP2A=3496490
k_ERBB4_PIK3=7
KM_ERBB4_PIK3=400000
k_cKIT_PIK3=7
KM_cKIT_PIK3=400000
k_PIK3_Ras=0.0771067
KM_PIK3_Ras=272056
k_Raf_Ras=0.884096
KM_Raf_Ras=62464.6
k_Mek_Raf=185.759
KM_Mek_Raf=4768350
k_Raf_Act=15.1212
KM_Raf_Act=119355
k_Ras_RasGab=1509.36
KM_Ras_RasGab=1432410
k_MEK_PP2A=2.83243
KM_MEK_PP2A=518753
k_Raf_RafPP=0.126329
KM_Raf_RafPP=1061.71
kd_P90Rsk=0.005
k_ERK_LKB1 =5.0
KM_ERK_LKB1=1500000.0
k_LKB1_AMPK=5
KM_LKB1_AMPK=1500000
k_AMPK_mTOR=7
KM_AMPK_mTOR=400000
k_Akt_mTOR=7
KM_Akt_mTOR=400000
k_mTOR_p70S6K=7
KM_mTOR_p70S6K=400000
k_Akt_AMPK=7
KM_Akt_AMPK=400000
k_Akt_Bad=7
KM_Akt_Bad=400000
k_Akt_CASP=7
KM_Akt_CASP=400000


#total proteins
Ras_tot=120000.0
Raf_tot=120000.0
Mek_tot=600000.0
Erk_tot=600000.0
P90Rsk_tot=120000.0
PIK3_tot=120000.0
Akt_tot=120000.0
LKB1_tot=360000.0
AMPK_tot=360000.0
mTOR_tot=360000.0
p70S6K_tot=360000.0
RafPP=120000
PP2A=120000
RasGapActive=120000


#initial conditions
EGFR_active_0=8000.0
IGFR_active_0=8000.0
ERBB4_active_0=8000
cKIT_active_0=8000
SOS_0=0.0
DSOS_0=120000.0
Ras_active_0=0.0
Raf_active_0=0.0
Mek_active_0=0.0
Erk_active_0=0.0
P90Rsk_active_0=0.0
PIK3_active_0=0.0
Akt_active_0=0.0
LKB1_active_0=0
AMPK_active_0=0
mTOR_active_0=0
p70S6K_active_0=0
BAD_active_0=120000
CASP_active_0=120000
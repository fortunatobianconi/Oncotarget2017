def case_report_model(x,t,p,Xt):

#kinetic parameters
    gamma_EGFR=p[0]
    gamma_IGF1R=p[1]
    gamma_EBB4=p[2]
    gamma_cKIT=p[3]
    kf_PI3K_active=p[4]
    k_p90Rsk_Erk=p[5]
    KM_p90Rsk_Erk=p[6]
    k_SOS_E=p[7]
    KM_SOS_E=p[8]
    k_SOS_ERBB4=p[9]
    KM_SOS_ERBB4=p[10]
    k_SOS_cKIT=p[11]
    KM_SOS_cKIT=p[12]
    k_Ras_SOS=p[13]
    KM_Ras_SOS=p[14]
    k_ERK_MEK=p[15]
    KM_ERK_MEK=p[16]
    k_DSOS_p90Rsk=p[17]
    KM_DSOS_p90Rsk=p[18]
    k_SOS_I=p[19]
    KM_SOS_I=p[20]
    k_PIK3_IGF1R=p[21]
    KM_PIK3_IGF1R=p[22]
    k_PIK3_EGFR=p[23]
    KM_PIK3_EGFR=p[24]
    k_Akt_PIK3=p[25]
    KM_Akt_PIK3=p[26]
    kd_Akt=p[27]
    k_Erk_PP2A=p[28]
    KM_Erk_PP2A=p[29]
    k_ERBB4_PIK3=p[30]
    KM_ERBB4_PIK3=p[31]
    k_cKIT_PIK3=p[32]
    KM_cKIT_PIK3=p[33]
    k_PIK3_Ras=p[34]
    KM_PIK3_Ras=p[35]
    k_Raf_Ras=p[36]
    KM_Raf_Ras=p[37]
    k_Mek_Raf=p[38]
    KM_Mek_Raf=p[39]
    k_Raf_Act=p[40]
    KM_Raf_Act=p[41]
    k_Ras_RasGab=p[42]
    KM_Ras_RasGab=p[43]
    k_MEK_PP2A=p[44]
    KM_MEK_PP2A=p[45]
    k_Raf_RafPP=p[46]
    KM_Raf_RafPP=p[47]
    kd_P90Rsk=p[48]
    k_ERK_LKB1=p[49]
    KM_ERK_LKB1=p[50]
    k_LKB1_AMPK=p[51]
    KM_LKB1_AMPK=p[52]
    k_AMPK_mTOR=p[53]
    KM_AMPK_mTOR=p[54]
    k_Akt_mTOR=p[55]
    KM_Akt_mTOR=p[56]
    k_mTOR_p70S6K=p[57]
    KM_mTOR_p70S6K=p[58]
    k_Akt_AMPK=p[59]
    KM_Akt_AMPK=p[60]
    k_Akt_Bad=p[61]
    KM_Akt_Bad=p[62]
    k_Akt_CASP=p[63]
    KM_Akt_CASP=p[64]
 
#total proteins
    Ras_tot=Xt[0]
    Raf_tot=Xt[1]
    Mek_tot=Xt[2]
    Erk_tot=Xt[3]
    P90Rsk_tot=Xt[4]
    PIK3_tot=Xt[5]
    Akt_tot=Xt[6]
    LKB1_tot=Xt[7]
    AMPK_tot=Xt[8]
    mTOR_tot=Xt[9]
    p70S6K_tot=Xt[10]
    RafPP=Xt[11]
    PP2A=Xt[12]
    RasGapActive=Xt[13]

#initial conditions
    EGFR_active=x[0]
    IGFR_active=x[1]
    ERBB4_active=x[2]
    cKIT_active=x[3]
    SOS=x[4]
    DSOS=x[5]
    Ras_active=x[6]
    Raf_active=x[7]
    Mek_active=x[8]
    Erk_active=x[9]
    P90Rsk_active=x[10]
    PIK3_active=x[11]
    Akt_active=x[12]
    LKB1_active=x[13]
    AMPK_active=x[14]
    mTOR_active=x[15]
    p70S6K_active=x[16]
    BAD_active=x[17]
    CASP_active=x[18]

#ODEs
    dEGFR_active=-gamma_EGFR*EGFR_active
    dIGFR_active=-gamma_IGF1R*IGFR_active
    dERBB4_active=-gamma_EBB4*ERBB4_active
    dcKIT_active=-gamma_cKIT*cKIT_active
    dSOS=k_SOS_ERBB4*ERBB4_active*DSOS/(KM_SOS_ERBB4+DSOS)+k_SOS_cKIT*cKIT_active*DSOS/(KM_SOS_cKIT+DSOS)+k_SOS_E*EGFR_active*DSOS/(KM_SOS_E+DSOS)+IGFR_active*k_SOS_I*DSOS/(KM_SOS_I+DSOS)-P90Rsk_active*k_DSOS_p90Rsk*SOS/(KM_DSOS_p90Rsk+SOS)
    dDSOS=-k_SOS_ERBB4*ERBB4_active*DSOS/(KM_SOS_ERBB4+DSOS)-k_SOS_cKIT*cKIT_active*DSOS/(KM_SOS_cKIT+DSOS)-k_SOS_E*EGFR_active*DSOS/(KM_SOS_E+DSOS)-IGFR_active*k_SOS_I*DSOS/(KM_SOS_I+DSOS) + P90Rsk_active*k_DSOS_p90Rsk*SOS/(KM_DSOS_p90Rsk+SOS)
    dRas_active=SOS*k_Ras_SOS*(Ras_tot-Ras_active)/(KM_Ras_SOS+(Ras_tot-Ras_active))-RasGapActive*k_Ras_RasGab*Ras_active/(KM_Ras_RasGab+Ras_active)
    dRaf_active=Ras_active*k_Raf_Ras*(Raf_tot-Raf_active)/(KM_Raf_Ras+(Raf_tot-Raf_active))-RafPP*k_Raf_RafPP*Raf_active/(KM_Raf_RafPP+Raf_active)-Akt_active*k_Raf_Act*Raf_active/(KM_Raf_Act+Raf_active)
    dMek_active=Raf_active*k_Mek_Raf*(Mek_tot-Mek_active)/(KM_Mek_Raf+(Mek_tot-Mek_active))-PP2A*k_MEK_PP2A*Mek_active/(KM_MEK_PP2A+Mek_active)
    dErk_active=Mek_active*k_ERK_MEK*(Erk_tot-Erk_active)/(KM_ERK_MEK+(Erk_tot-Erk_active))-PP2A*k_Erk_PP2A*Erk_active/(KM_Erk_PP2A+Erk_active)
    dP90Rsk_active=Erk_active*k_p90Rsk_Erk*(P90Rsk_tot-P90Rsk_active)/(KM_p90Rsk_Erk+(P90Rsk_tot-P90Rsk_active))-kd_P90Rsk*P90Rsk_active
    dPIK3_active=IGFR_active*k_PIK3_IGF1R*(PIK3_tot-PIK3_active)/(KM_PIK3_IGF1R+(PIK3_tot-PIK3_active))+k_PIK3_EGFR*EGFR_active*(PIK3_tot-PIK3_active)/(KM_PIK3_EGFR+(PIK3_tot-PIK3_active))+k_ERBB4_PIK3*ERBB4_active*(PIK3_tot-PIK3_active)/(KM_ERBB4_PIK3+(PIK3_tot-PIK3_active))+k_cKIT_PIK3*cKIT_active*(PIK3_tot-PIK3_active)/(KM_cKIT_PIK3+(PIK3_tot-PIK3_active))+Ras_active*k_PIK3_Ras*(PIK3_tot-PIK3_active)/(KM_PIK3_Ras+(PIK3_tot-PIK3_active))-kf_PI3K_active*PIK3_active
    dAkt_active=PIK3_active*k_Akt_PIK3*(Akt_tot-Akt_active)/(KM_Akt_PIK3+(Akt_tot-Akt_active))-kd_Akt*Akt_active
    dLKB1_active=k_ERK_LKB1*Erk_active*(LKB1_tot-LKB1_active)/(KM_ERK_LKB1+(LKB1_tot-LKB1_active))
    dAMPK_active=k_LKB1_AMPK*(LKB1_tot-LKB1_active)*(AMPK_tot-AMPK_active)/(KM_LKB1_AMPK+(AMPK_tot-AMPK_active))-k_Akt_AMPK*Akt_active*AMPK_active/(KM_Akt_AMPK+AMPK_active)
    dmTOR_active=k_Akt_mTOR*Akt_active*(mTOR_tot-mTOR_active)/(KM_Akt_mTOR+(mTOR_tot-mTOR_active))-k_AMPK_mTOR*AMPK_active*mTOR_active/(KM_AMPK_mTOR+mTOR_active)
    dp70S6K_active=k_mTOR_p70S6K*mTOR_active*(p70S6K_tot-p70S6K_active)/(KM_mTOR_p70S6K+(p70S6K_tot-p70S6K_active))
    dBAD_active=-k_Akt_Bad*Akt_active*BAD_active/(KM_Akt_Bad+BAD_active)
    dCASP_active=-k_Akt_CASP*Akt_active*CASP_active/(KM_Akt_CASP+CASP_active)
    

    
    return [dEGFR_active, dIGFR_active, dERBB4_active, dcKIT_active, dSOS, dDSOS, dRas_active, dRaf_active, dMek_active, dErk_active, dP90Rsk_active, dPIK3_active, dAkt_active, dLKB1_active, dAMPK_active, dmTOR_active, dp70S6K_active, dBAD_active, dCASP_active]

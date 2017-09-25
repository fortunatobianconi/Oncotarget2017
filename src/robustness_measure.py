import numpy as np
import pyDOE
from memmap import memmap
from joblib import Parallel,delayed
import parameters as parameters
from model_sim import model_sim
import pickle

  
#initial conditions
x0=[parameters.EGFR_active_0, parameters.IGFR_active_0, parameters.ERBB4_active_0, parameters.cKIT_active_0, parameters.SOS_0,
     parameters.DSOS_0, parameters.Ras_active_0, parameters.Raf_active_0, parameters.Mek_active_0, parameters.Erk_active_0,
     parameters.P90Rsk_active_0, parameters.PIK3_active_0, parameters.Akt_active_0, parameters.LKB1_active_0,
     parameters.AMPK_active_0, parameters.mTOR_active_0, parameters.p70S6K_active_0, parameters.BAD_active_0,
     parameters.CASP_active_0]
         

#kinetic parameters
p=[parameters.gamma_EGFR, parameters.gamma_IGF1R, parameters.gamma_ERBB4, parameters.gamma_cKIT, 
   parameters.kf_PI3K_active, parameters.k_p90Rsk_Erk, parameters.KM_p90Rsk_Erk, parameters.k_SOS_E, 
   parameters.KM_SOS_E, parameters.k_SOS_ERBB4, parameters.KM_SOS_ERBB4, parameters.k_SOS_cKIT, parameters.KM_SOS_cKIT,
   parameters.k_Ras_SOS, parameters.KM_Ras_SOS, parameters.k_ERK_MEK, parameters.KM_ERK_MEK, parameters.k_DSOS_p90Rsk,
   parameters.KM_DSOS_p90Rsk, parameters.k_SOS_I, parameters.KM_SOS_I, parameters.k_PIK3_IGF1R, parameters.KM_PIK3_IGF1R, 
   parameters.k_PIK3_EGFR, parameters.KM_PIK3_EGFR, parameters.k_Akt_PIK3, parameters.KM_Akt_PIK3, parameters.kd_Akt, 
   parameters.k_Erk_PP2A, parameters.KM_Erk_PP2A, parameters.k_ERBB4_PIK3, parameters.KM_ERBB4_PIK3, parameters.k_cKIT_PIK3,
   parameters.KM_cKIT_PIK3, parameters.k_PIK3_Ras, parameters.KM_PIK3_Ras, parameters.k_Raf_Ras, parameters.KM_Raf_Ras, 
   parameters.k_Mek_Raf, parameters.KM_Mek_Raf, parameters.k_Raf_Act, parameters.KM_Raf_Act, parameters.k_Ras_RasGab,
   parameters.KM_Ras_RasGab, parameters.k_MEK_PP2A, parameters.KM_MEK_PP2A, parameters.k_Raf_RafPP, parameters.KM_Raf_RafPP,
   parameters.kd_P90Rsk, parameters.k_ERK_LKB1, parameters.KM_ERK_LKB1, parameters.k_LKB1_AMPK, parameters.KM_LKB1_AMPK, 
   parameters.k_AMPK_mTOR, parameters.KM_AMPK_mTOR, parameters.k_Akt_mTOR, parameters.KM_Akt_mTOR, parameters.k_mTOR_p70S6K,
   parameters.KM_mTOR_p70S6K, parameters.k_Akt_AMPK, parameters.KM_Akt_AMPK, parameters.k_Akt_Bad, parameters.KM_Akt_Bad,
   parameters.k_Akt_CASP, parameters.KM_Akt_CASP]

#total proteins
Xt=[parameters.Ras_tot, parameters.Raf_tot, parameters.Mek_tot, parameters.Erk_tot, parameters.P90Rsk_tot, 
    parameters.PIK3_tot, parameters.Akt_tot, parameters.LKB1_tot, parameters.AMPK_tot, parameters.mTOR_tot,
    parameters.p70S6K_tot, parameters.RafPP, parameters.PP2A, parameters.RasGapActive]

#lower bound
LBpi=0.1

#upper bound
UBpi=10

#number of realizations
Nr=100

#number of samples of the parameter vector
NSample=100000

#definition of the number of proteins to measure
ProteinNumber=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
Results_tot=np.zeros(shape=(Nr,len(ProteinNumber),NSample))

#parameters not to be perturbed
fixed_p=[]
fixed_Xt=[]

#initial conditions to perturb
unfixed_x0=[0, 1, 2, 3, 5, 17, 18]

Np=len(p)+len(Xt)+len(unfixed_x0)-len(fixed_p)-len(fixed_Xt)
PerturbationT=np.zeros(shape=(Nr,NSample,Np))
PerturbationPT=np.zeros(shape=(Nr,NSample,len(p)-len(fixed_p)))
PerturbationXtT=np.zeros(shape=(Nr,NSample,len(Xt)-len(fixed_Xt)))
PerturbationXT=np.zeros(shape=(Nr,NSample,len(unfixed_x0)))

for k in range(0,Nr):
    print k

    #perturbation of the parameter vector using latin hypercube sampling
    PerturbationP=LBpi+(UBpi-LBpi)*pyDOE.lhs((np.asarray(p).size-len(fixed_p)),samples=NSample) 
    PerturbationX=LBpi+(UBpi-LBpi)*pyDOE.lhs(len(unfixed_x0),samples=NSample)
    PerturbationXt=LBpi+(UBpi-LBpi)*pyDOE.lhs((np.asarray(Xt).size-len(fixed_Xt)),samples=NSample)
    Perturbation=np.concatenate(([PerturbationP],[PerturbationX],[PerturbationXt]),axis=2)


    Results=np.zeros((len(ProteinNumber),NSample))
       
    all_data=memmap([ProteinNumber,p,PerturbationP,fixed_p, Xt,PerturbationXt,fixed_Xt, x0,PerturbationX, unfixed_x0])

    #parallel model simulation, njobs is the number of workers
    Output=Parallel(n_jobs=10,verbose=9,batch_size='auto',max_nbytes='1K')(delayed(model_sim)(all_data[0],all_data[1],all_data[2],all_data[3], all_data[4],all_data[5],all_data[6],all_data[7],all_data[8],all_data[9], i) for i in range(NSample))

    for j in range(0,len(ProteinNumber)):
        Results[j,:]=([Output[v][j] for v in range(0,NSample)])

    Results_tot[k,:,:]=Results
    PerturbationPT[k,:,:]=PerturbationP
    PerturbationXtT[k,:,:]=PerturbationXt
    PerturbationXT[k,:,:]=PerturbationX
    PerturbationT[k,:,:]=Perturbation
    Perturbation=np.concatenate(([PerturbationP],[PerturbationX],[PerturbationXt]),axis=2)

#save results
with open('FileName.pickle', 'w') as f:
    pickle.dump([Results_tot, PerturbationT, PerturbationPT, PerturbationXtT, PerturbationXT, p, x0, Xt, fixed_p,fixed_Xt,unfixed_x0,Nr, NSample], f)
    
   
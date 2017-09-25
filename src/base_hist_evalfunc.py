from __future__ import division
from upper_lower_set import upper_lower_set
from intersection import intersection
import pickle
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from collections import defaultdict
import os

#set current directory
path=os.getcwd()
os.chdir(path)
total_path= path + '/FileName.pickle'

#load measure
with open(total_path) as f:
    Results_tot, PerturbationT, PerturbationPT, PerturbationXtT, PerturbationXT, p, x0, Xt, fixed_p, fixed_Xt, unfixed_x0, Nr, NSample = pickle.load(f)

# add 1 where parameters are fixed
for i in range(0, len(fixed_p)):
    PerturbationPT = np.insert(PerturbationPT, fixed_p[i], 1, axis=2)

for i in range(0, len(fixed_Xt)):
    PerturbationXtT = np.insert(PerturbationXtT, fixed_Xt[i], 1, axis=2)

PerturbationT = np.concatenate((PerturbationPT, PerturbationXT, PerturbationXtT), axis=2)


#definition of the names of proteins to evaluate
protein_name=['Erk','AMPK','mTOR']

#threshold definition for upper and lower tail
lowThr=[]
highThr=[]

#upper and lower tail construction
indexes=upper_lower_set(Results_tot,Nr,protein_name,lowThr,highThr)

#intersection among proteins
region_tot=defaultdict(list)
for k in range(0,Nr):
    [region, region_name]=intersection(indexes[k])

    for i in range(0,len(region_name)):    
        region_tot[region_name[i]].append(region[i])



MIRIT=[]
CloudH_T=[]
CloudL_T=[]

#for each realization
for k in range(0,Nr):
    
    #cloud definition
    NCloudH=0
    NCloudL=0
    
    #specify the target and opposite behavior for the selected proteins for calibration
    NCloudH=len(region_tot['ULUL'][k])
    NCloudL=len(region_tot['LULU'][k])
    NCloud=min([NCloudH, NCloudL])

    X=[x0[0],x0[1],x0[3]]
    pN=p+X+Xt     
    
    MIRI=[]
    CloudH=[]
    CloudL=[]

    #for each parameter evaluate the probability density function
    for ip in range(0,len(pN)-1):
        print ip

        #skip parameters whose value is fixed
        if(PerturbationT[k,0,ip]!=1):
            start=np.amin(pN[ip]*PerturbationT[k,:,ip])
            stop=np.amax(pN[ip]*PerturbationT[k,:,ip])
            step=pN[ip]*(NCloud**(-1./3))
            parameter_axis=np.linspace(start,stop,(stop-start)/step)
            parameter_axis = (parameter_axis[:, np.newaxis])

        #realizations of the parameter vector giving the desired region
            H=pN[ip]*PerturbationT[k,region_tot['ULL'][k],ip]
            L=pN[ip]*PerturbationT[k,region_tot['LUU'][k],ip]
       
            x_H=H[:,np.newaxis]
    
            std_H=np.std(x_H)
            bw_H=(4*(std_H)**5/(3*len(x_H)))**(0.2)  #bandwidth estimation
    
            kde_H=KernelDensity(kernel='gaussian',bandwidth=bw_H).fit(x_H)
            log_H=kde_H.score_samples(parameter_axis) #conditioned probability density function estimation

            x_L=L[:,np.newaxis]
    
            std_L=np.std(x_L)
            bw_L=(4*(std_L)**5/(3*len(x_L)))**(0.2)  #bandwidth estimation
    
            kde_L=KernelDensity(kernel='gaussian',bandwidth=bw_L).fit(x_L)
            log_L=kde_L.score_samples(parameter_axis) #conditioned probability density function estimation

            kde_H_list=np.exp(log_H).tolist()
            kde_L_list=np.exp(log_L).tolist()
    
            kde_H_I=kde_H_list.index(max(np.exp(log_H)))
            kde_L_I=kde_L_list.index(max(np.exp(log_L)))

            #plot pdf
            plt.figure(ip)
            plt.plot(parameter_axis,np.exp(log_H))
            plt.plot(parameter_axis,np.exp(log_L))
            plt.show()

        #save parameter values
            CloudH.append(parameter_axis[kde_H_I])
            CloudL.append(parameter_axis[kde_L_I])
     
        #calculate MIRI
            array_H=np.array(kde_H_list)/sum(kde_H_list)
            array_L=np.array(kde_L_list)/sum(kde_L_list)
            s_Xi_ks = sum(abs(array_H - array_L))
            MIRI.append(s_Xi_ks)
     
        
        CloudH_T.append(CloudH)
        CloudL_T.append(CloudL)
        MIRIT.append(MIRI)

#boxplot MIRI for all realizations
plt.figure()
plt.boxplot(np.asarray(MIRIT))
plt.show()

#bar MIRI of the last realization
x_axis=np.linspace(1,len(MIRI),len(MIRI))
plt.figure()
plt.bar(x_axis,MIRI)
plt.show()  
     

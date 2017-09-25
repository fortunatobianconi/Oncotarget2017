# This function computes the evaluation function of each protein in ProteinNumber for each sample i
# of the parameter vector

import numpy as np
from scipy.integrate import odeint
from case_report_model import case_report_model


def model_sim(ProteinNumber, p, PerturbationP, fixed_p, Xt, PerturbationXt, fixed_Xt, x0, PerturbationX, unfixed_x0, i):

    # add 1 where parameters are fixed
    PerturbationP=PerturbationP.tolist()
    for j in range(0,len(PerturbationP)):  
        for s in range(0,len(fixed_p)):
            PerturbationP[j].insert(fixed_p[s],1)
    pi=p*np.asarray(PerturbationP[i][:])
      
    PerturbationXt=PerturbationXt.tolist()
    for j in range(0,len(PerturbationXt)):  
        for s in range(0,len(fixed_Xt)):
            PerturbationXt[j].insert(fixed_p[s],1)
    Xti=Xt*np.asarray(PerturbationXt[i][:])
    
    y_0i=list(x0)
    t=0
    for j in unfixed_x0:
        y_0i[j]=PerturbationX[i][t]*x0[j]
        t=t+1

    #time definition
    t =np.linspace(0,2000,125000)

    #model integration
    y = odeint(case_report_model, y_0i, t, args=(pi, Xti), mxstep=500000)

    #measure of the evaluation function chosen
    areas=[]
    for t in range(0,len(ProteinNumber)):
        area=0
        for j in range(0,len(y)):
            area=area+y[j][ProteinNumber[t]]
        areas.append(area)

    return areas

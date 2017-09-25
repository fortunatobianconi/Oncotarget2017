# for each protein in input, it returns the set of realizations of the parameter vector for which
# the protein evaluation function has either an high or low value. It takes in input:
# - Results: the set of evaluation functions measured
# - Nr: number of realizations
# - NSample: number of samples of the parameter vector
# - proteinName: name of measured proteins
# -lowThr, highThr: specify how to select the high or low tail of a protein. It is possible to
# choose a different threshold for each protein


def upper_lower_set(Results,Nr,proteinName,lowThr,highThr):
    
    #set to nan possible negative values in the measured results

    indexset_tot = list()

    for k in range(0,Nr):
        for j in range(0,len(Results[k])):
            for i in range(0,len(Results[k,j])):
                if(Results[k][j][i] < 0):
                    Results[k][j][i] = float('nan')

    #for each realization
    for k in range(0,Nr):
        index_set=list()
        for s in range(0,len(proteinName)):  #for each protein

            #definition of the arrays that contain the indexes for the upper and lower tail
            upper_name = 't_'+proteinName[s]+'_Upper'
            lower_name = 't_'+proteinName[s]+'_Lower'
            mydict={upper_name:[],lower_name:[]}

            for t in range(0,len(Results[k][s])):  #scan measures of each protein
                if Results[k][s][t] > highThr[s]:
                    mydict[upper_name].append(t)
                elif Results[k][s][t] < lowThr[s]:
                    mydict[lower_name].append(t)   

            #insert 0 if there are no indexes for the desired region
            if mydict[upper_name] == []:
                mydict[upper_name] = list([0])
                      
            if mydict[lower_name] == []:
                mydict[lower_name] = list([0])
            
            #the two sets are bound in a single array
            index_set.append(mydict[upper_name])
            index_set.append(mydict[lower_name])

        indexset_tot.append(index_set)

    return indexset_tot
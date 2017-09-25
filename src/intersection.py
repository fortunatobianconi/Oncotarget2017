# intersection is a recursive function. It takes in input the array of indexes returned by
# the upper_lower_set function and intersects each single set with the others. In this way, it finds
# all the samples of the parameter vector that give simultaneously the desired behavior for all proteins measured.


def intersection(indexes):

    region = list()
    region_name = list()

    #base case
    if(len(indexes) == 4):
        region.append(list(set(indexes[0]).intersection(set(indexes[2]))))
        region_name.append('UU')
               
        region.append(list(set(indexes[0]).intersection(set(indexes[3]))))
        region_name.append('UL')
       
        region.append(list(set(indexes[1]).intersection(set(indexes[2]))))
        region_name.append('LU')
       
        region.append(list(set(indexes[1]).intersection(set(indexes[3]))))
        region_name.append('LL')
     
        return region, region_name

    #recursive case
    else:
        
        partial_region=intersection(indexes[0:len(indexes)-2])
     
        for i in range (0,len(partial_region[0])):
            region.append(list(set(partial_region[0][i]).intersection(set(indexes[len(indexes)-2]))))
            region_name.append(partial_region[1][i]+'U')

            region.append(list(set(partial_region[0][i]).intersection(set(indexes[len(indexes)-1]))))
            region_name.append(partial_region[1][i]+'L')

        return region, region_name
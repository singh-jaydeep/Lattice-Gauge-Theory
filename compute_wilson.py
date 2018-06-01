import constants as c
import numpy as np
import array_manip as am
import gen_purpose_funs as gpf

## We only do positive orientation loops, so only need the desired line

def compute_wilson(location_array,time_dist,sp_dist,line,field_array):

    location1 = location_array
    #print(location1)
    loop = 1.0

    for i in range(0,time_dist):
        loop = gpf.gen_mult_2(loop, am.array_retrieve(location1,0,1,field_array))
        location1 = am.move_location(location1,0,1,field_array) ## move timelike
    for i in range(0,sp_dist):
        loop = gpf.gen_mult_2(loop, am.array_retrieve(location1,line,1,field_array))
        location1 = am.move_location(location1,line,1,field_array)
    for i in range(0,time_dist):
        loop = gpf.gen_mult_2(loop, am.array_retrieve(location1,0,0,field_array))
        location1 = am.move_location(location1,0,0,field_array)
    for i in range(0,sp_dist):
        loop = gpf.gen_mult_2(loop, am.array_retrieve(location1,line,0,field_array))
        location1 = am.move_location(location1,line,0,field_array)

    if c.gaugegp == 1:
        return loop
    else:
        return np.trace(loop)

def global_wilson_meas(field_array):
    array = np.zeros((c.ran_time-1,c.ran_spat-1),dtype=np.complex_)
    shape = np.shape(field_array)


    if c.dim == 2:
        shape = shape[0:2]
    if c.dim == 3:
        shape = shape[0:3]
    if c.dim == 4:
        shape = shape[0:4]



    for index in np.ndindex((c.ran_time-1,c.ran_spat-1)):
        counter = 0.0
        temp = 0.0
        line = np.random.randint(1,c.dim)
        for location in np.ndindex(shape):
            temp = temp + compute_wilson(location,index[0]+1,index[1]+1,line,field_array)
            counter =  counter + 1.0
            #for line in range(1,c.dim):  ## one side is already in time direction, so need others to be diff
            #    temp = temp + compute_wilson(location,index[0]+1,index[1]+1,line,field_array)
            #counter = counter + (c.dim-1)
        array[index[0],index[1]] = temp/counter

    return np.reshape(array,(1,(c.ran_time-1)*(c.ran_spat-1)))


## Computes the mean value of the plaquettes in a given field configuration
def compute_plaq_mean(field_array):
    shape = np.shape(field_array)

    ## extract location information
    if c.dim == 2:
        shape = shape[0:2]
    if c.dim == 3:
        shape = shape[0:3]
    if c.dim == 4:
        shape = shape[0:4]

    temp = 0
    counter = 0
    #for location in np.nindex(shape):

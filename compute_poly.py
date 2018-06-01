import array_manip as am
import constants as c
import gen_purpose_funs as gpf
import numpy as np


## For now, we are arbitrarily picking a non-time direction of the field array (always
## the second direction), and computing the correlator between points at location
## [0,0] and [0,k] as k ranges between 0 and the lattice size/2 (b/c periodic boundary conditions
## render latticesize/2 as the largest distance between two loops)
### IF IN NEED OF MORE DATA POINTS, ADD TO THIS FUNCTION
def global_corr_measure(field_array):

    if c.dim==2:
        location1 = [0,0]
        location2 = [0,0]
        array = []
        temp = []

        for i in range(0,c.latt_size):
            for j in range(0,int(c.latt_size/2)+1):
                temp = np.append(temp,correlator(location1,location2,field_array))
                location2 = am.move_location(location2,1,1,field_array)
            if i == 0:
                array = np.copy(temp)
            else:
                array = np.vstack((array,temp))
            temp = []
            location1 = am.move_location(location1,1,1,field_array)
            location2 = am.move_location(location1,1,1,field_array)

    if c.dim ==3:
        location1 = [0,0,0]
        location2 = [0,0,0]
        array = []
        temp = []

        for index in np.ndindex((1,c.latt_size,c.latt_size)):
            location1 = index
            location2 = index
            for l in range(0,int(c.latt_size/2)+1):
                temp = np.append(temp,correlator(location1,location2,field_array))
                location2 = am.move_location(location2,1,1,field_array)
            if len(array) == 0:
                array = np.copy(temp)
            else:
                array = np.vstack((array,temp))
            temp = []

    elif c.dim == 4:
        location1 = [0,0,0,0]
        location2 = [0,0,0,0]
        array = []
        temp = []

        for index in np.ndindex((1,c.latt_size,c.latt_size,c.latt_size)):
            location1 = index
            location2 = index
            for l in range(0,int(c.latt_size/2)+1):
                temp = np.append(temp,correlator(location1,location2,field_array))
                location2 = am.move_location(location2,1,1,field_array)
            if len(array) == 0:
                array = np.copy(temp)
            else:
                array = np.vstack((array,temp))
            temp = []

    return array


def poly_loop(location_array,field_array):
    #print(location_array)
    current_loc = np.copy(location_array)
    line = 0 ## the time dimension is implicitly represented as the first dimension
    dirc = 1 ## we always move in the same direction to compute the Polyakov loop
    loop_prod = 1

    for i in range(0,c.time_latt_size):
        loop_prod = gpf.gen_mult_2(loop_prod,am.array_retrieve(current_loc,line,dirc,field_array))
        current_loc = am.move_location(current_loc,line,dirc,field_array)

    if c.gaugegp == 1:
        return loop_prod
    else:
        return np.trace(loop_prod)

def correlator(location1,location2,field_array):
    return poly_loop(location1,field_array)*gpf.herm_conj(poly_loop(location2,field_array))

import numpy as np
import constants as c
import gen_purpose_funs as gpf
import pickle
import os
import compute_poly as poly


def initialize():
    ##idea: Put everything in a massive array of dimensions
    ##(latt_size)^dim. At each lattice point, have an array of size
    ##(2)^dim/2, with jth entry the matrix for the jth direction from
    ##that point. Storing things in all directions (so added redundancy)

    ## U(1) easy because each link element is just a complex number
    if c.gaugegp == 1 :
        if c.dim == 2:
            field_array = np.ones((c.time_latt_size,c.latt_size,c.dim,2,1),dtype=np.complex_)
        if c.dim == 3:
            field_array = np.ones((c.time_latt_size,c.latt_size,c.latt_size,c.dim,2,1),dtype=np.complex_)
        if c.dim == 4:
            field_array = np.ones((c.time_latt_size,c.latt_size,c.latt_size,c.latt_size,c.dim,2,1),dtype=np.complex_)
            print(np.shape(field_array))
        return field_array
    elif c.gaugegp == 2:
        if c.dim == 2:
            field_array = np.ones((c.time_latt_size,c.latt_size,c.dim,2,2,2),dtype=np.complex_)
        if c.dim ==3:
            field_array = np.ones((c.time_latt_size,c.latt_size,c.latt_size\
                          ,c.dim,2,2,2),dtype=np.complex_)
        if c.dim == 4:
            field_array = np.ones((c.time_latt_size,c.latt_size,c.latt_size,c.latt_size\
                          ,c.dim,2,2,2),dtype=np.complex_)

        ## Need to initialize the field as identity matrices
        for index in np.ndindex(np.shape(field_array)):
            if index[c.dim+2] == index[c.dim+3]: ## these are the diagonals
                field_array[index] = 1.0
            else:
                field_array[index] = 0.0

        return field_array
    elif c.gaugegp == 3:
        if c.dim == 2:
            field_array = np.ones((c.time_latt_size,c.latt_size,c.dim,2,3,3),dtype=np.complex_)
        if c.dim == 3:
            field_array = np.ones((c.time_latt_size,c.latt_size,c.latt_size\
                          ,c.dim,2,3,3),dtype=np.complex_)
        if c.dim == 4:
            field_array = np.ones((c.time_latt_size,c.latt_size,c.latt_size,c.latt_size\
                          ,c.dim,2,3,3),dtype=np.complex_)

        ## Need to initialize the field as identity matrices
        for index in np.ndindex(np.shape(field_array)):
            if index[c.dim+2] == index[c.dim+3]: ## these are the diagonals
                field_array[index] = 1.0
            else:
                field_array[index] = 0.0

        return field_array


## A centralized way to extract elements of the field_array without needing to specify the dimension
## Implicitly, we use the convention that a 1 in dirc corresponds to "up/right", 0 to "left/down"
def array_retrieve(location_array,line,dirc,field_array):

    if c.gaugegp ==1:
        if c.dim == 2:
            return field_array[location_array[0],location_array[1],line,dirc,0]
        if c.dim == 3:
            return field_array[location_array[0],location_array[1],location_array[2],line,dirc,0]
        elif c.dim ==4:
            return field_array[location_array[0],location_array[1],location_array[2],location_array[3],line,dirc,0]
    if c.gaugegp == 2 or c.gaugegp == 3:
        if c.dim == 2:
            return field_array[location_array[0],location_array[1],line,dirc,:,:]
        if c.dim == 3:
            return field_array[location_array[0],location_array[1],location_array[2],\
            line,dirc,:,:]
        if c.dim == 4:
            return field_array[location_array[0],location_array[1],location_array[2],\
            location_array[3],line,dirc,:,:]


## Given a point in the lattice and a direction, returns the adjacent point, incorporating
## periodic boundary conditions. line refers to the dimension along which we are moving
## and dir = +-1 indicates if the entry is increasing or decreasing
def move_location(location_array,line,dirc,field_array):
    new_location = np.copy(location_array)

    #print('location moving from is with line', location_array, line,dirc)
    ## deal with time direction boundary condition
    if line == 0 and new_location[0] == (c.time_latt_size-1) and dirc == 1:
        new_location[0] = 0
        return new_location
    elif line == 0 and new_location[0] == 0 and  dirc == 0:
        new_location[0] = c.time_latt_size-1
        return new_location

    ## deal with spatial direction boundary condition
    if new_location[line] == (c.latt_size-1) and dirc == 1:
        new_location[line] = 0
    elif new_location[line] == 0 and dirc == 0:
        new_location[line] = c.latt_size-1
    ## if not at any boundary, just add
    else:
        if dirc == 1:
            new_location[line] = new_location[line] + 1
        if dirc == 0:
            new_location[line] = new_location[line] - 1
    return new_location

## Given the location of the link you want to update, returns a field
## with the updated link, adjusting the adjacent opposite orientation link
## so that the condition b(x,y) = b(y,x)^-1 is always met
def array_change_link(location_array,line,dirc,field_array,newlink):
    #temp_array = np.copy(field_array)
    temp_array = field_array
    temp_location = np.copy(location_array)

    if c.gaugegp == 1:
        inverselink = 1/newlink
        if c.dim == 2:
            temp_array[location_array[0],location_array[1],line,dirc,0] = newlink
            nextlocation = move_location(temp_location,line,dirc,field_array)
            temp_array[nextlocation[0],nextlocation[1],line,1-dirc,0] = inverselink
        if c.dim == 3:
            temp_array[location_array[0],location_array[1],location_array[2], \
            line,dirc,0] = newlink
            nextlocation = move_location(temp_location,line,dirc,field_array)
            temp_array[nextlocation[0],nextlocation[1],location_array[2], \
                line,1-dirc,0] = inverselink
        elif c.dim == 4:
            temp_array[location_array[0],location_array[1],location_array[2], \
            location_array[3], line,dirc,0] = newlink
            nextlocation = move_location(temp_location,line,dirc,field_array)
            temp_array[nextlocation[0],nextlocation[1],location_array[2], \
                location_array[3],line,1-dirc,0] = inverselink
    elif c.gaugegp == 2 or c.gaugegp == 3:
        inverselink = gpf.herm_conj(newlink)
        if c.dim == 2:
            temp_array[location_array[0],location_array[1],line,dirc,:,:] = newlink
            nextlocation = move_location(temp_location,line,dirc,field_array)
            temp_array[nextlocation[0],nextlocation[1],line,1-dirc,:,:] = inverselink
        if c.dim == 3:
            temp_array[location_array[0],location_array[1],location_array[2],\
            line,dirc,:,:] = newlink
            nextlocation = move_location(temp_location,line,dirc,field_array)
            temp_array[nextlocation[0],nextlocation[1],nextlocation[2],\
            line,1-dirc,:,:] = inverselink
        if c.dim == 4:
            temp_array[location_array[0],location_array[1],location_array[2],\
            location_array[3],line,dirc,:,:] = newlink
            nextlocation = move_location(temp_location,line,dirc,field_array)
            temp_array[nextlocation[0],nextlocation[1],nextlocation[2],\
            nextlocation[3],line,1-dirc,:,:] = inverselink

    return temp_array


def load():
    (run_string,data_string) = gpf.gen_string()

    if os.path.isfile(run_string):
        with open(run_string, 'rb') as pickle_in:
            field_array = pickle.load(pickle_in)
            print("Loading equilibriated array")
    else:
        print("No array present, generating new array")
        field_array = initialize()

    if os.path.isfile(data_string):
        with open(data_string, 'rb') as pickle_in:
            correlator_array = pickle.load(pickle_in)
            print("Loading correlator array")
    else:
        print("No correlator array present, generating new array")
        correlator_array = []#poly.global_corr_measure(field_array)

    return (field_array,correlator_array)

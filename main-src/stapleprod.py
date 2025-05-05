import constants as c
import numpy as np
import array_manip as am
import gen_purpose_funs as gpf


##Compute a staple for a link in a given new direction, along with the reflected
##staple with opposite orientation. The total staple is the sum of all indstaples
##across new directions not equal to the direction passed in

## To define a staple, we need two lines (defining a plane)
def indstaple(location_array,line,line_perp,direction,field_array):
    ## Initialize things

    if c.gaugegp == 1:
        stapleproduct1 = 1  # positive orientation staple
        stapleproduct2 = 1    # negative orientation staple
    elif c.gaugegp == 2:
        stapleproduct1 = np.array([[1,0],[0,1]])

        stapleproduct2 = np.array([[1,0],[0,1]])
    elif c.gaugegp == 3:
        stapleproduct1 = np.array([[1,0,0],[0,1,0],[0,0,1]])
        stapleproduct2 = np.array([[1,0,0],[0,1,0],[0,0,1]])


    ## Maybe I can be super explicit with this, given the work done to make field_array as
    ## transparent as possible (5/4/18)

    newposition1_1 = am.move_location(location_array,line,direction,field_array)
    newposition2_1= am.move_location(newposition1_1,line_perp,direction,field_array)
    newposition3_1 = am.move_location(newposition2_1,line,1-direction,field_array)

    #print(stapleproduct1,am.array_retrieve(newposition1_1,line_perp,direction,field_array))
    stapleproduct1 = gpf.gen_mult_4(stapleproduct1 \
        ,am.array_retrieve(newposition1_1,line_perp,direction,field_array)\
        ,am.array_retrieve(newposition2_1,line,1-direction,field_array)\
        ,am.array_retrieve(newposition3_1,line_perp,1-direction,field_array))


    newposition1_2 = am.move_location(location_array,line,direction,field_array)
    newposition2_2 = am.move_location(newposition1_2,line_perp,1-direction,field_array)
    newposition3_2 = am.move_location(newposition2_2,line,1-direction,field_array)


    stapleproduct2 = gpf.gen_mult_4(stapleproduct2 \
        ,am.array_retrieve(newposition1_2,line_perp,1-direction,field_array)\
        ,am.array_retrieve(newposition2_2,line,1-direction,field_array)\
        ,am.array_retrieve(newposition3_2,line_perp,direction,field_array))


    return stapleproduct1 + stapleproduct2

def stapleproduct(location_array,line,direction,field_array):

    sum = 0
    for i in range(0,c.dim):
        if i != line:
            sum = sum + indstaple(location_array,line,i,direction,field_array)

    return sum

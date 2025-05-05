## This function accepts a configuration, and the random bags input, randomly selects one lattice node
## and direction, updating the link via a random bag element, and accepting according to the Metropolis
## selection rule.

## For the selection, given the candidate link, call a function to compute the staple product for
## this node/direction. Pass the staple product, and new/old configuration to a function to compute
## delta S, returning a boolean (accept/don't accept). Repeat this repeat_iters_node times, before
## returning. Make sure Iter counter is aware that more than 1 iter has passed during a call to
## update_config

##


import numpy as np
import constants as c
import random
import stapleprod
import array_manip as am
import math
import compute_poly as poly
import gen_purpose_funs as gpf


def sweep_update(field_array,randombag,algorithm):
    accept_counter = 0
    total_counter = 0
    shape = np.shape(field_array)

    if c.gaugegp == 1:
        index_shape = shape[0:(len(shape)-1)]
    else:
        index_shape = shape[0:(len(shape)-2)]

    for index in np.ndindex(index_shape):
        location = index[0:c.dim]

        #print(location)
        line = index[c.dim]
        direction = index[c.dim+1]
        link = am.array_retrieve(location,line,direction,field_array)
        #print(link)
        staple = stapleprod.stapleproduct(location,line,direction,field_array)
        #print(location)

        if algorithm == "Metropolis":
            for j in range(0,c.repeatnode_iters):
                temp = random.choice(randombag)
                #print(temp)
                #print(link)
                #print(gpf.gen_mult_2(temp,link))
                newlink = gpf.make_unitary(gpf.gen_mult_2(temp,link))
                #print('newlink is ', newlink)
                total_counter += 1
                if accept_update(staple,link, newlink):
                    accept_counter += 1
                    #print('accepted update')
                    field_array = am.array_change_link(location,line,direction,field_array,newlink)


        if algorithm == "HeatBath":

            ## Generates random numbers for heat bath update
            a = np.sqrt(np.linalg.det(staple))
            for j in range(0,c.repeatnode_iters):
                total_counter +=1
                accept_counter += 1 ## Note we are accepting everything, so the accept rate isn't super relevant
                accepted1 = 0
                l_squared = 0
                while accepted1 == 0:
                    randnum = np.random.rand(1,3)
                    randnum = 1 - randnum[0]
                    l_squared = -1 / (2*c.coupling*a) * (np.log(randnum[0])+math.cos(2*np.pi*randnum[1])**2 \
                    *np.log(randnum[2]))
                    if (np.random.rand())**2 <= 1-l_squared:
                        accepted1 = 1
                x0 = 1 - 2*l_squared
                magx = np.sqrt(1-(x0)**2)

                accepted2 = 0
                x1 = 0
                x2 = 0
                x3 = 0
                while accepted2 == 0:
                    randnum = np.random.rand(1,3)
                    randnum = 2*randnum[0] - 1 ## three numbers in [-1,1]
                    if np.linalg.norm(randnum) <= 1:
                        accepted2 = 1
                        magr = np.sqrt(randnum[0]**2+randnum[1]**2+randnum[2]**2)
                        x1 = randnum[0]/magr * magx
                        x2 = randnum[1]/magr * magx
                        x3 = randnum[2]/magr * magx

                temp = gpf.gen_SU2_mat(x0,x1,x2,x3)
                newlink = gpf.make_unitary(gpf.gen_mult_2(temp,gpf.herm_conj(staple/a)))
                #newlink = gpf.make_unitary(temp,link)
                #temp = gpf.make_unitary(gpf.gen_mult_2(gpf.gen_SU2_mat(x0,x1,x2,x3),link))
                #newlink = gpf.gen_mult_2(temp,gpf.herm_conj(staple/a))
                field_array = am.array_change_link(location,line,direction,field_array,newlink)
        #break
    return (field_array,accept_counter,total_counter)

def accept_update(staple,link,newlink):
    #Computes change in action
    if c.gaugegp == 1:
        tr_del_S = -1.0*c.coupling*np.real((newlink-link)*staple)
    elif c.gaugegp == 2:
        #print("staple is", staple)
        #print("link is", link)
        #print("newlink is", newlink)
        #print("newlink-link is", newlink-link)
        tr_del_S = -1/2.0*c.coupling*np.real(np.trace(gpf.gen_mult_2((newlink-link),staple)))
    elif c.gaugegp == 3:
        tr_del_S = -1/3.0*c.coupling*np.real(np.trace(gpf.gen_mult_2((newlink-link),staple)))
    ## We always accept if we are strictly decreasing the action (tr_del_S<0)
    ## but probabalistically accept configurations increasing action. This is
    ## what allows us to explore the space of field configurations

    if np.random.rand() < math.exp(-1*tr_del_S):
        return 1
    else:
        return 0

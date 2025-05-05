## Builds a random array of matrices in a given gauge group
import constants as c
import gen_purpose_funs as gpf
import numpy as np
import random

def randombag(group,randomint):  ## We pass in the group, because the SU(3) case needs to call the SU(2) case
    if group == 1:
        u1_array = np.random.rand(1,c.size_rand_bag)
        u1_array = u1_array[0]
        u1_array = 2*randomint*u1_array - randomint
        #u1_array = (1 + 1j*u1_array) / np.sqrt(1+u1_array*u1_array)
        u1_array = np.exp(1j*u1_array)
        u1_conj = np.conj(u1_array)
         ## to ensure symmetric selection probability,
        ## the random matrix must have the inverse of every group element
        return np.append(u1_array,u1_conj)

    elif group == 2:
        ## for each random bag element, start with four numbers distributed in [-1/2,1/2]
        ## for proof of this, see page 83 of Gattringer and Lang
        su2_bag = []
        collection1 = np.random.rand(3,c.size_rand_bag) - 1/2 ## spatial components


        collection1 = randomint*collection1/np.linalg.norm(collection1,axis=0)

        collection2 = np.random.rand(1,c.size_rand_bag) - 1/2 ## time component

        collection2 = collection2[0]
        collection2 = np.sign(collection2)*np.sqrt(1-(randomint)**2)


        for elem in range(0,c.size_rand_bag):
            temp1 = gpf.gen_SU2_mat(collection2[elem],collection1[0,elem],collection1[1,elem],collection1[2,elem])
            temp2 = np.conj(np.transpose(temp1))
            su2_bag.append(temp1)
            su2_bag.append(temp2)
        return su2_bag

    elif group == 3:
        su3_bag = []
        bag1 = np.random.permutation(randombag(2,randomint))
        bag2 = np.random.permutation(randombag(2,randomint))
        bag3 = np.random.permutation(randombag(2,randomint))  ## Generate three bags of SU(2) random variables



        for i in range(0,len(bag1)):

            ## Generate embedded matrix1
            mat1 = np.pad(bag1[i], (0,1), 'constant', constant_values=0)
            mat1[2,2] = 1

            ## Generate embedded matrix 2
            mat2 = np.pad(bag2[i], (1,0), 'constant', constant_values=0)
            mat2[0,0] = 1

            ## Generate embedded matrix3
            mat3 = np.zeros((3,3),dtype=np.complex_)
            mat3[0,0] = bag3[i][0,0]
            mat3[0,2] = bag3[i][0,1]
            mat3[2,0] = bag3[i][1,0]
            mat3[2,2] = bag3[i][1,1]
            mat3[1,1] = 1

            ## Generate random SU(3) matrix and its inverse
            mat_su3 = gpf.gen_mult_3(mat1,mat3,mat2)
            mat_su3_inv = np.conj(np.transpose(mat_su3))
            su3_bag.append(mat_su3)
            su3_bag.append(mat_su3_inv)


        return su3_bag

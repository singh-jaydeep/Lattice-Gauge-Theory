import numpy as np
import constants as c

## Multiplies both numbers and matrices in a unified manner
def gen_mult_2(a,b):
    return np.dot(a,b)

def gen_mult_3(a,b,c):
    return np.dot(a,np.dot(b,c))

def gen_mult_4(a,b,c,d):
    return np.dot(a, np.dot(b, np.dot(c,d)))

## Not really necessary,tbh
def herm_conj(a):
    if c.gaugegp ==1:
        return np.conj(a)
    else:
        return np.conj(np.transpose(a))

## Given a, outputs a unitary matrix according to a prescription
def make_unitary(a):
    if c.gaugegp == 1:
        return a/np.absolute(a)
    elif c.gaugegp == 2:
        row1 = a[[0],:]/np.linalg.norm(a[[0],:])
        row1 = row1[0]
        row2 = np.array([-1*np.conj(row1[1]),np.conj(row1[0])])
        comp = np.array([row1,row2])
        return comp
    elif c.gaugegp == 3:
        row1 = a[[0],:]/np.linalg.norm(a[[0],:])
        row1 = row1[0]
        row2 = a[[1],:] - row1*np.dot(a[[1],:],np.conj(row1))
        row2 = row2 / np.linalg.norm(row2)
        row2 = row2[0]
        row3 = np.cross(np.conj(row1),np.conj(row2))
        comp = np.array([row1,row2,row3])
        return comp

def gen_string():
    string1 = f"field_arrays/equil_gp{c.gaugegp}_dim{c.dim}_coup{c.coupling}_latt{c.latt_size}_time{c.time_latt_size}_alg{c.algorithm}_obs{c.observable}"
    string2 = f"correlation_arrays/equil_gp{c.gaugegp}_dim{c.dim}_coup{c.coupling}_latt{c.latt_size}_time{c.time_latt_size}_alg{c.algorithm}_obs{c.observable}_DATA"
    return (string1,string2)

def group_to_string():
    if c.gaugegp == 1:
        return "U(1)"
    elif c.gaugegp == 2:
        return "SU(2)"
    elif c.gaugegp == 3:
        return "SU(3)"

## Uses the pauli matrices to generate a SU(2) matrix given parameters
## Four parameters together should be norm 1
def gen_SU2_mat(x0,x1,x2,x3):
    pauli1 = np.array([[0,1],[1,0]])
    pauli2 = np.array([[0,-1*1j],[1j,0]])
    pauli3 = np.array([[1,0],[0,-1]])
    ident = np.array([[1,0],[0,1]])

    return x0*ident + 1j * (x1*pauli1) + 1j * (x2*pauli2) + 1j * (x3*pauli3)

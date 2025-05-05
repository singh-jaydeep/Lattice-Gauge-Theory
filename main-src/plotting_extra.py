import constants as c
import randombag as rm
import array_manip as am
import numpy as np
import output
import gen_purpose_funs as gpf
import exec_simulation_run as run
import matplotlib.pyplot as plt
import stapleprod as sp
import update_config as update
import pickle
import random
import os
import analyze_data as analyze
import compute_poly as poly
import compute_wilson as wilson
import time


pickle_in = open("correlation_arrays/equil_gp1_dim2_coup4_latt16_time16_algMetropolis_obsWilson_DATA","rb")
corr_matrix1 = pickle.load(pickle_in)

pickle_in = open("correlation_arrays/equil_gp1_dim2_coup5_latt16_time16_algMetropolis_obsWilson_DATA","rb")
corr_matrix2 = pickle.load(pickle_in)

pickle_in = open("correlation_arrays/equil_gp1_dim2_coup6_latt16_time16_algMetropolis_obsWilson_DATA","rb")
corr_matrix3 = pickle.load(pickle_in)

(a1,b1,c1) = analyze.process_updated(corr_matrix1,5000)
(a2,b2,c2) = analyze.process_updated(corr_matrix2,5000)
(a3,b3,c3) = analyze.process_updated(corr_matrix3,5000)

plt.errorbar(a1,b1,yerr = c1, label ="beta = 4")
plt.errorbar(a2,b2,yerr = c2,label ="beta = 5")
plt.errorbar(a3,b3,yerr = c3,label ="beta = 6")
plt.legend()
plt.xlabel('Particle Distance a|m|')
plt.ylabel('aV(a|m|)')
plt.title('Static Quark/Anti-Quark Potential U(1)')
plt.show()

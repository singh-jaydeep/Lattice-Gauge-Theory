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

pickle_in = open("correlation_arrays/equil_gp3_dim2_coup4_latt16_time16_algMetropolis_obsWilson_DATA","rb")
corr_matrix= pickle.load(pickle_in)
shape = np.shape(corr_matrix)
print(shape)
shaved_corr = corr_matrix[(shape[0]-300):shape[0],:]
analyze.process(np.ones((1,1)),corr_matrix,1,1,1)
#analyze.process_updated(shaved_corr,5000)
sample_points = []
for i in range(0,len(corr_matrix)):
     temp = corr_matrix[i]
     sample_points.append(temp[24])
print(np.abs(sample_points))
analyze.draw_traceplot(corr_matrix,5,5)

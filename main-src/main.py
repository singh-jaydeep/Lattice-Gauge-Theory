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



(A,B,C,D) = run.simulation_run_sweeps()



#####################################################################################
######################## Sample data analysis code ##################################
#####################################################################################
#pickle_in = open("correlation_arrays/equil_gp1_dim4_coup0.98_latt16_time16_algMetropolis_obsWilson_DATA","rb")
#corr_matrix= pickle.load(pickle_in)
#shape = np.shape(corr_matrix)
#shaved_corr = corr_matrix[(shape[0]-60):shape[0],:]
#analyze.process(np.ones((1,1)),shaved_corr,1,1,1)

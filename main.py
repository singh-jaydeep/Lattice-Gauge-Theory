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

#field = am.initialize()
#print(field)
#print(wilson.global_wilson_meas(field))
#print(wilson.compute_wilson([0,0,0,0],4,4,1,field))
#print(np.shape(field))
#print("retrieved is", am.array_retrieve([0,0,0],1,1,field))
#bag = rm.randombag(3,c.random_spread)
#print(bag)
#print(field)
#test = random.choice(bag)
#print(test)
#print(np.dot(test,gpf.herm_conj(test)))
#(field,counter,counter2) = update.sweep_update(field,bag,"Metropolis")
#print('done with one iteration')
#print(np.shape(field))
#print(field[0,0,0,0,:,:])
#update.update_config(field,bag)
#wilson.global_wilson_meas(field)
#print(poly.poly_loop([0,0],field))
#print(poly.correlator([0,0],[0,1],field))

#A = np.array([[1,2,3],[2,3,4],[4,5,6]])
#B = gpf.make_unitary(A)
#print(np.dot(B,np.transpose(B)))

(A,B,C,D) = run.simulation_run_sweeps()
#(A,B,C,D) = run.simulation_run()



#####################################################################################
#pickle_in = open("correlation_arrays/equil_gp1_dim4_coup0.98_latt16_time16_algMetropolis_obsWilson_DATA","rb")
#corr_matrix= pickle.load(pickle_in)
#shape = np.shape(corr_matrix)
#print(shape)
#shaved_corr = corr_matrix[(shape[0]-60):shape[0],:]
#analyze.process(np.ones((1,1)),shaved_corr,1,1,1)
#analyze.process_updated(shaved_corr,5000)
#sample_points = []
#for i in range(0,len(corr_matrix)):
#     temp = corr_matrix[i]
#     sample_points.append(temp[2])
#print(np.abs(sample_points))
#temp1 = np.abs(np.reshape(corr_matrix[1,:],(c.ran_time-1,c.ran_spat-1)))
#temp2 = np.abs(np.reshape(corr_matrix[10,:],(c.ran_time-1,c.ran_spat-1)))
#np.reshape(corr_matrix[1,:],(c.ran_time-1,c.ran_spat-1))
#print(temp1[0:6,0:6])
#print(temp2[0:6,0:6])
#print(np.reshape(corr_matrix[1,:],(7,5)))
#print(np.shape(corr_matrix))
#pickle_in = open("field_arrays/equil_gp3_dim4_coup5.7_latt8_time6_algMetropolis","rb")
#field = pickle.load(pickle_in)
#wilson.global_wilson_meas(field)
#print("done")

#analyze.process(np.ones((1,1)),corr_matrix,1,1)
#corr_matrix = corr_matrix[0]
#print(corr_matrix[0][0])
#shape = np.shape(corr_matrix)
#indx = 30*np.arange(0,int(shape[0]/30))
#corr_matrix = corr_matrix[0:shape[0],1:shape[1]]
#corr_matrix = corr_matrix[indx,1:shape[1]]
#print(field[0:3,0:3,0,0,0])
#print(corr_matrix[0:10,0:4])
#analyze.bootstrap_analyze(corr_matrix,10)
#analyze.bootstrap_analyze(corr_matrix[(shape[0]-50000):shape[0],:],5000)











#print(np.shape(corr_array))

#print(np.shape(field))
#print(np.shape(poly.global_corr_measure(field)))
#print(poly.global_corr_measure(field))
#shape = np.shape(corr_matrix)
#corr_clean = corr_matrix[1:shape[0],:]
#analyze.process(np.zeros((1,1)),corr_clean,0,1)

#(A,B,C,D) = run.simulation_run()

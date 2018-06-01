import constants as c
import numpy as np
import randombag
import compute_poly as poly
import gen_purpose_funs as gpf
import array_manip as am
import update_config as update
import stapleprod as staple
import output
import analyze_data as analyze
import matplotlib.pyplot as plt
import compute_wilson as wilson
import pickle

def simulation_run_sweeps():
    random_scale = c.random_spread
    accepted_counter = 0
    totalupdate_counter = 0

    (field,correlator_array) = am.load() ## uses pre-equilibriated lattice if one exists, otherwise generates new array

    bag = randombag.randombag(c.gaugegp,random_scale)
    iters = 0 ## Note about iters: we are counting each iteration of the below loop, which in fact
              ## updates the configuration c.repeatnode_iters times. So iters is not the total number
              ## of updates to the field, but the number of updates modulo these repetitive updates
    meas_bool = 0 ## Are we recording measurements yet
    completed_sweeps = 0
    completed_meas = 0

    sweep_vec = []
    test_vec = []  ## For eventual trace plot

    if c.dim == 2:
        origin = [0,0]
    elif c.dim == 3:
        origin = [0,0,0]
    else:
        origin = [0,0,0,0]

    while iters < c.total_iters_sweeps:

        if meas_bool == 0:
            (field,accepted,total) = update.sweep_update(field,bag,c.algorithm)
            accepted_counter = accepted_counter + accepted
            totalupdate_counter = totalupdate_counter + total

            completed_sweeps += 1
            print("Completed equilibriation sweep", completed_sweeps)
            if completed_sweeps == c.equil_sweeps:
                meas_bool = 1

        #test_vec.append(wilson.compute_wilson(origin,3,3,1,field))
        #sweep_vec.append(iters)

        if iters % c.inter_rand_iters == 0:
            bag = randombag.randombag(c.gaugegp,random_scale)

        print("accept rate is currently", accepted_counter/totalupdate_counter)

        if meas_bool == 1 and c.observable == "Polyakov":
            if iters % c.inter_meas_sweeps == 0:
                new_meas = poly.global_corr_measure(field)
                completed_meas += 1
                print("Completed ", completed_meas, "/", c.total_measurements, "measurements")
                if len(correlator_array) == 0:
                    correlator_array = new_meas
                else:
                    correlator_array = np.vstack((correlator_array,new_meas))
            else:
                (field,accepted,total)= update.sweep_update(field,bag,c.algorithm)
                accepted_counter = accepted_counter + accepted
                totalupdate_counter = totalupdate_counter + total
                print("Completed intermediate sweep")

        if meas_bool == 1 and c.observable == "Wilson":
            if iters % c.inter_meas_sweeps == 0:
                new_meas =wilson.global_wilson_meas(field)
                completed_meas += 1
                print("Completed ", completed_meas, "/", c.total_measurements, "measurements")
                if len(correlator_array) == 0:
                    correlator_array = new_meas
                else:
                    correlator_array = np.vstack((correlator_array,new_meas))
            else:
                (field,accepted,total)= update.sweep_update(field,bag,c.algorithm)
                accepted_counter = accepted_counter + accepted
                totalupdate_counter = totalupdate_counter + total
                print("Completed intermediate sweep")

        if iters == 0 or iters % c.update_user_iters == 0:
            output.updateuser(iters,meas_bool)

        if iters % c.inter_save_iters == 0:
            print("Saving")
            (run_string,data_string) = gpf.gen_string()
            print("Updating field stored in ", run_string)
            with open(run_string, 'wb') as pickle_out:
                pickle.dump(field,pickle_out)
                pickle_out.close()

            print("Updating correlator stored in ", data_string)
            with open(data_string, 'wb') as pickle_out:
                pickle.dump(correlator_array,pickle_out)
                pickle_out.close()

        iters += 1

    (dist_vec,corr_analyzed) = analyze.process(field,correlator_array,accepted_counter,totalupdate_counter,random_scale)

    (run_string,data_string) = gpf.gen_string()
    print("Updating field stored in ", run_string)
    with open(run_string, 'wb') as pickle_out:
        pickle.dump(field,pickle_out)
        pickle_out.close()

    print("Updating correlator stored in ", data_string)
    with open(data_string, 'wb') as pickle_out:
        pickle.dump(correlator_array,pickle_out)
        pickle_out.close()

    #plt.plot(sweep_vec,test_vec)
    #plt.xlabel("Sweep Number")
    #plt.ylabel("Wilson Loop Value 3x3")
    #plt.title("MCMC Trace Plot " +"Inv Coupling = " + str(c.coupling) + " Dim = " + str(c.dim))
    #plt.show()

    return (field,correlator_array,dist_vec,corr_analyzed)

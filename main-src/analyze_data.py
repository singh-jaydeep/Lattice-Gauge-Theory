import constants as c
import matplotlib.pyplot as plt
import gen_purpose_funs as gpf
import numpy as np
from scipy.optimize import curve_fit
import pylab



## Underlying physics being used here: formula -log(<P(m)P(n)>) ~ C*latt_size*V(r)
def process(field_array,corr_matrix,accepted,totalupdates,randomscale):

    if c.observable == "Polyakov":
        shape = np.shape(corr_matrix)
        print(shape)
        corr_clean = corr_matrix[1:shape[0],:] ## take out the trivial first row
        corr_aver = np.abs(np.sum(corr_clean,axis=0)) / (shape[0]-1)
        corr_log_scale = -1*np.log(corr_aver) /c.time_latt_size
        #print(corr_log_scale[1])
        #print(corr_log_scale[2])
        #print(corr_log_scale[3])

        dist_vec = np.arange(0,shape[1])

        #### TEST PART OF code
        #### Prints a graph of the variance in a data point with respect to num of measurements

        test_var = []
        #for i in range(1,shape[0]-1):
        #    temp = np.std(np.abs(corr_clean[1:i+1,3]),axis=0)*1/np.sqrt(i) / (np.abs(np.sum(corr_clean[1:i+1,3],axis=0)) / (i))
        #    test_var.append(temp)
        #plt.plot(np.arange(0,len(test_var)),test_var)
        #plt.ylabel('variance vs measurement number')
        #plt.show()


        #######

        ## Computes error bars, assuming uncorrelated data
        corr_dev = np.std(np.abs(corr_clean),axis=0)

        result = corr_log_scale

        print(np.max(1/np.sqrt(shape[0])*corr_dev))
        print(np.max(corr_dev))
        print(np.min(corr_aver))
        print(np.max(corr_aver))
        print(np.min(corr_matrix))
        print(np.max(corr_matrix))
        print('the accept rate of this run was', accepted/totalupdates)
        print('the starting random spread statistic was ', c.random_spread, "and ended on ", randomscale)

        group = gpf.group_to_string()
        plt.errorbar(dist_vec,corr_log_scale, yerr=1/np.sqrt(shape[0]-1)*corr_dev/corr_aver)
        plt.xlabel('Distance between loops a|m|')
        plt.ylabel('aV(a|m|)')
        plt.title('Static Quark/Anti-Quark Potential' + group)
        plt.show()

    ## Handed in a
    if c.observable == "Wilson":
        print('the accept rate of this run was', accepted/totalupdates)
        ## unpackage the rows of the matrix:
        shape_orig = np.shape(corr_matrix)  ## e.g. N * 48
        print(shape_orig)
        running_total = np.zeros((c.ran_time-1,c.ran_spat-1))
        counter = 0

        for i in range(0,shape_orig[0]):
            running_total = running_total + np.reshape(corr_matrix[i,:],(c.ran_time-1,c.ran_spat-1))
            counter = counter+1
            #corr_revised.append(np.reshape(corr_matrix[i,:],(c.ran_time-1,c.ran_spat-1)))
                    ## the minus ones are appearing b/c we don't want to include loops of zero area. So
                    ## the 0th element corresponds to one unit in that direction.


        av_sample = np.abs(running_total/counter)
        #av_sample = np.abs(np.real(running_total/counter))


        ###### But our samples also have lots of errors, so let's use loops up to 6x6
        #av_sample = av_sample[0:6,0:6]
        #print(np.log(av_sample))
        new_shape = np.shape(av_sample)
        print(new_shape)
        time_ran = new_shape[0]
        spat_ran = new_shape[1]



        potential = []
        t_vec = 1+np.arange(time_ran)

        ## At fixed r, how does it depend on t?

        potential_test = []
        for i in range(0,time_ran):
            m = -1*np.log(av_sample[i,3])
            potential_test.append(m)

        plt.plot(t_vec,potential_test)
        plt.title("-ln(W(R,T)) vs T with fixed R")
        plt.show()
        ##########################################

        for i in range(0,spat_ran): ## For each spat_distance r, we will compute V(r) for that distance
            #print(np.log(av_sample[2:8,i]))
            #m,b = pylab.polyfit(t_vec, -1*np.log(av_sample[:,i]), 1)
            #m,b = pylab.polyfit(t_vec[0:3], -1*np.log(av_sample[0:3,i]), 1)

            #(param,corr) = curve_fit(functional_form,dist_vec,potential)


            m =  -1*np.log(av_sample[6,i])/6
            #m = -1*np.log(av_sample[5,i]) - -1*np.log(av_sample[3,i])
            potential.append(m)

        result = potential
        dist_vec = np.arange(1,spat_ran+1)
        #dist_vec = np.arange(1,5)+1

        group = gpf.group_to_string()
        plt.plot(dist_vec,potential)
        plt.xlabel('Distance between loops a|m|')
        plt.ylabel('aV(a|m|)')
        plt.title('Static Quark/Anti-Quark Potential' + group)
        plt.show()

    return (dist_vec,result)







def bootstrap_analyze(corr_matrix,boot_iters):
    shape = np.shape(corr_matrix)
    corr_matrix_new = corr_matrix[1:shape[0],:]

    new_sample = []

    for i in range(0,boot_iters):
        idx = np.random.randint(0, shape[0]-1, (1,shape[0]-1))
        idx = idx[0]
        copy = np.copy(corr_matrix[idx,:])
        corr_aver = np.abs(np.sum(copy,axis=0)) / (shape[0]-1)
        corr_log_scale = -1*np.log(corr_aver) /c.time_latt_size
        if i == 0:
            new_sample = corr_log_scale
        else:
            new_sample = np.vstack((new_sample,corr_log_scale))

    ##Computes new mean
    sample_av = np.mean(new_sample,axis=0)
    ##Computes new error
    sample_err = ( 1/(boot_iters)*np.sum( (new_sample-sample_av)**2,axis=0) )**(.5)

    group = gpf.group_to_string()
    dist_vec = np.arange(0,shape[1])
    scale = compute_scale(dist_vec,sample_av)
    plt.errorbar(dist_vec,sample_av,yerr = sample_err)
    plt.axvline(x=.5/scale,linestyle = '--')
    plt.xlabel('Distance between loops a|m|')
    plt.ylabel('aV(a|m|)')
    plt.title('Static Quark/Anti-Quark Potential' + group)
    plt.legend(["Potential",".5 fm Sommer Parameter"])
    plt.show()

    return (dist_vec,sample_av,sample_err)

def functional_form(n,A,B,C):
    return A + B/n + C*n

## takes the potential V(r) as input, outputs the length scale of a in femtometers(10^-15)
def compute_scale(dist_vec,potential):

    (param,corr) = curve_fit(functional_form,dist_vec,potential)
    Bopt = param[1]
    Copt = param[2]
    scale = .5 * np.sqrt((Copt / (1.65+Bopt)))
    print("One lattice unit corresponds to a scale of ", scale, "femtometers")
    #dist_vec_physical = scale * dist_vec
    #plt.errorbar(dist_vec_physical,sample_av,yerr = sample_err)
    return scale



def process_updated(corr_matrix,boot_iters):
    shape_orig = np.shape(corr_matrix)
    new_samples = []

    for i in range(0,boot_iters):
        #running_total = np.zeros((c.ran_time-1,c.ran_spat-1))
        running_total_first = np.zeros((1,shape_orig[1]))
        counter = 0

        idx = np.random.randint(0, shape_orig[0]-1, (1,shape_orig[0]))
        idx = idx[0]
        copy = np.copy(corr_matrix[idx,:])
        for i in range(0,shape_orig[0]):
            running_total_first = running_total_first + copy[i,:]
            counter = counter+1
        av_sample = np.abs(running_total_first/counter)
        print(np.shape(av_sample),np.shape(new_samples))
        if len(new_samples) == 0:
            new_samples = av_sample
        else:
            new_samples = np.vstack((new_samples,av_sample))


    # From the above samples, we need to compute V(R) from each sample
    len_range = np.arange(1,c.ran_spat+1)
    time_vec = np.arange(1,c.ran_time)
    results_mat = np.zeros((boot_iters,len(len_range)))
    for i in range(0,c.ran_spat-1):
        for j in range(0,boot_iters):
            temp = np.reshape(new_samples[j,:],(c.ran_time-1,c.ran_spat-1))
            #m,b = pylab.polyfit(time_vec, -1*np.log(temp[:,i]), 1)
            #results_mat[j,i] = m
            results_mat[j,i] = -np.log(temp[4,i])/4

    group = gpf.group_to_string()
    mean_potential = np.mean(results_mat,axis=0)
    error_bars = np.std(results_mat,axis=0)
    plt.errorbar(len_range[0:7],mean_potential[0:7],yerr = error_bars[0:7])
    plt.xlabel('Distance between loops a|m|')
    plt.ylabel('aV(a|m|)')
    plt.title('Static Quark/Anti-Quark Potential ' + group)
    plt.show()

    return(len_range[0:7],mean_potential[0:7],error_bars[0:7])

def draw_traceplot(corr_matrix,time_size,spat_size): ## which loop are we tracking
    tracked = []
    shape = np.shape(corr_matrix)
    num_samples = shape[0]
    for i in range(0,num_samples):
        reshaped = np.reshape(corr_matrix[i,:],(c.ran_time-1,c.ran_spat-1))
        tracked.append(np.abs(reshaped[time_size-1,spat_size-1]))

    plt.plot(np.arange(0,num_samples),tracked)
    plt.xlabel("Measurement Number")
    plt.ylabel("Wilson Loop Value")
    plt.title("MCMC Trace Plot " +"Inv Coupling = " + str(c.coupling) + " Dim = " + str(c.dim))
    plt.show()

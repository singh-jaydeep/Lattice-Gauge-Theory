#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#### Extra code from update_config.py after altering field_array to separately
#### include the line (dimension) under consideration, and the +-1 direction

location = np.random.randint(0,c.latt_size-1,size=(1,4)) ## need up to 4 numbers
direction = np.random.randomint(0,2*c.dim)

  ## Case if dim == 2
  if c.dim == 2:
      ## Identify link under consideration
      link = field_array(location[0],location[1],direction,0)
      ## If working with U(1)

      if c.gaugegp == 1:
          ## For speed, we can compute the staple prior to iteration
          ## Write this!!
          staple = stapleprod.stapleproduct([location[0],location[1]],direction,field_array)
          for j in range(0,c.repeat_iters_node):

              newlink = random.choice(randombag)*link
              ## For U(1), projecting to unitarity is easy
              newlink = newlink/abs(newlink)
              ## Write the accept function!
              if accept_update(staple,link, newlink):
                  field_array(location[0],location[1],direction,0) = newlink
                  ### Update the opposite orientation link from the next point!!

      elif c.gaugegp == 2:
          ## do nothing yet

      return

  ## Case if dim == 4
  if c.dim == 4:
      return field_array

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################


##### Test code for variance
pickle_in = open("correlation_arrays/equil_gp1_dim2_coup0.5_latt20_time10_DATA","rb")
corr_matrix = pickle.load(pickle_in)

shape = np.shape(corr_matrix)
corr_clean = np.real(corr_matrix[1:shape[0],:])
#corr_clean = np.imag(corr_matrix[1:shape[0],:])
 ## take out the trivial first row
corr_aver = np.abs(np.sum(corr_clean,axis=0)) / (shape[0]-1)
#corr_aver = corr_aver[1:len(corr_aver)]
print(corr_aver)
corr_log_scale = -1*np.log(corr_aver) #/c.latt_size


dist_vec = np.arange(0,shape[1])

#### TEST PART OF code
#### Prints a graph of the variance in a data point with respect to num of measurements

test_var = []
for i in range(shape[0]-10000,shape[0]-1):
    temp = np.std(np.abs(corr_clean[1:i+1,3]),axis=0)*1/np.sqrt(i) / (np.abs(np.sum(corr_clean[1:i+1,3],axis=0)) / (i))
    test_var.append(temp)
plt.plot(np.arange(0,len(test_var)),test_var)
plt.ylabel('variance vs measurement number')
plt.show()



################################################################
#More variance test code
#####################################################################################
pickle_in = open("correlation_arrays/equil_gp1_dim2_coup2_latt20_time10_DATA","rb")
corr_matrix = pickle.load(pickle_in)

shape = np.shape(corr_matrix)
corr_clean = np.real(corr_matrix[1:shape[0],:])
#corr_clean = np.imag(corr_matrix[1:shape[0],:])
 ## take out the trivial first row
corr_aver = np.abs(np.sum(corr_clean,axis=0)) / (shape[0]-1)
#corr_aver = corr_aver[1:len(corr_aver)]
print("corr av is", corr_aver)
corr_log_scale = -1*np.log(corr_aver) #/c.latt_size

corr_dev = np.std(np.abs(corr_clean),axis=0)
print("corr dev is", corr_dev)
yerr=1/np.sqrt(shape[0]-1)*corr_dev/corr_aver
print("square root n factor is",1/np.sqrt(shape[0]-1) )
print("yerror is", yerr)
########################################################

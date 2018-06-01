coupling = 4
gaugegp = 3## 1 = U(1), 2 = SU(2), 3 = SU(3)
dim = 2 ## dimension of space-time
latt_size = 16## extent of spatial lattice in dimensionless lattice units
time_latt_size = 16## 10 is a test number in general

latt_points = (latt_size)**(dim-1) * time_latt_size
num_bonds = latt_points*(2*dim) ## total number of links
num_dirs = 2*dim ## from a given lattice point, the number of directions

size_rand_bag=200

## Controls the spread (since SU(3) randombag computed using SU(2) bags, adjust su2_randint for gauge gps 2,3)
#u1_randint=.5
#su2_randint = .06
random_spread = .05  ### RETURN TO 3.0



## Update algorithm
algorithm = "Metropolis" ## Options include "Metropolis" (For all gauge groups) and "HeatBath" (for SU(2) currently)

## What observable are you measuring?
observable = "Wilson" ## Options inclde "Wilson" or "Polyakov"

## Wilson specific (SWITCH BACK TO 10,10)
ran_time = 9
ran_spat = 9

###Constants that mark important iteration numbers/periods
repeatnode_iters= 1 ## Number of times that a particular node/direction is updated in a row
inter_rand_iters = 50 ## Should be <1/2 of the size of the random_bag. Number of iters between calls to refresh
                      ## the random bag
update_user_iters = 1000  ## How often to update the user
inter_save_iters = 10

total_measurements = 200 #300
equil_sweeps = 1#50
inter_meas_sweeps = 10
check_accept_rate_sweeps = 2
total_iters_sweeps = equil_sweeps + total_measurements*inter_meas_sweeps

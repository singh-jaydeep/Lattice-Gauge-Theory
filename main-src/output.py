import constants as c

def updateuser(iters, meas_bool):
    if iters == 0:
        print('Beginning simulation run with gauge group', c.gaugegp,\
            'and coupling', c.coupling)
    elif meas_bool:
        print('On iteration', iters, 'recording measurements')
    else:
        print('On iteration', iters, 'NOT recording measurements')

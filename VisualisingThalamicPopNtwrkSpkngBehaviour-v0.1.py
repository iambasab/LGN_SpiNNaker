''' 
LOG 27TH SEP BASAB: CODE MODIFIED BY AKASH TO COMPUTE MEAN OF MEAN MEMBRANE POTENTIALS GENERATED ACROSS MULTIPLE SIMULATIONS.

LOG 15th September Basab: Using the code that Akash wrote last year to visualise the data generated with SysManCyberB_Toni_Design_September16.py - looking at the correlation of the spike raster of the three LGN cell populations.

Code Copyright of Akash Bhattacharya, 2nd Year Student, MSc Course in Chemistry with Molecular Physics, Imperial College London, 2015-2016.
'''

#!/usr/bin/python
import numpy as np
import matplotlib.pylab as plt
from pylab import *
from pyNN.random import NumpyRNG, RandomDistribution


#def check_if_len_zero(*args):
#    verdict = [True] * len(args)
#    for count in range(0, len(args)):
#        if len(args[count]) > 0:
#            verdict[count] = False
#        else:
#            print("WARNING! The input index number %i has no results!", count)
#    return verdict


# Declaring constants.
TimeInt = 10            # milliseconds
TotalDuration = 2000   # milliseconds; total simulation time
TotalDataLength = TotalDuration * TimeInt


scale_fact = 10
NumCellsTCR = 5 * scale_fact    # Number of Cells of TCR
NumCellsTRN = 4 * scale_fact    # Number of Cells of TRN
NumCellsIN = 1 * scale_fact     # Number of Cells in IN


loop = input('Enter the number of loops in your simulation: \n')

''' Creating zeroes-array to store mean values'''
TCR_mean = np.array([[0] * TotalDataLength] * loop)
TRN_mean = np.array([[0] * TotalDataLength] * loop)
IN_mean = np.array([[0] * TotalDataLength] * loop)

# Repeating for all files.


for j in range(0, loop):
    '''Load Data'''
    # Each simulation of the SysManCyberB_Tonic_Design_September16.py file is run for 5 times.
    # So Loop can be any number from 1 to 5

    

    #Loading files
    n = 1  # input('enter number of simulation')
    f = 8   # input('enter frequency of simulation')
    foldername="Sim%d_%dhz_0916" % (n, f)

#    spike_source_ex = np.loadtxt('./'+foldername+'/spikesource_ex_'+`loop`+'.dat')
#    spike_source_inh = np.loadtxt('./'+foldername+'/spikesource_inh_'+`loop`+'.dat')
#    TCRspikes = np.loadtxt('./'+foldername+'/TCRspikes_'+`loop`+'.dat')
#    INspikes = np.loadtxt('./'+foldername+'/INspikes_'+`loop`+'.dat')
#    TRNspikes = np.loadtxt('./'+foldername+'/TRNspikes_'+`loop`+'.dat')


    TCRmempot = np.loadtxt('./'+foldername+'/TCRmempot_'+`loop`+'.dat')
    INmempot = np.loadtxt('./'+foldername+'/INmempot_'+`loop`+'.dat')
    TRNmempot = np.loadtxt('./'+foldername+'/TRNmempot_'+`loop`+'.dat')

    # Sanity check: do results exist or not?
#    check_results = check_if_len_zero(spike_source_ex, spike_source_inh, TCRspikes, INspikes, TRNspikes)
    

    '''MATRICES FORMATION'''
    signalTCR = []
    signalTRN = []
    signalIN = []

    # Creating lists of relevant data.
    for i in range(0, len(TCRmempot)):
        signalTCR.append(TCRmempot[i][2])
    for i in range(0, len(TRNmempot)):
        signalTRN.append(TRNmempot[i][2])
    for i in range(0, len(INmempot)):
        signalIN.append(INmempot[i][2])

    # Turning lists into arrays and reshaping into matrices.
    signalTCR = np.array(signalTCR)
    signalTRN = np.array(signalTRN)
    signalIN = np.array(signalIN)

    TCR = signalTCR.reshape(NumCellsTCR, TotalDataLength)
    IN = signalIN.reshape(NumCellsIN, TotalDataLength)
    TRN = signalTRN.reshape(NumCellsTRN, TotalDataLength)

    '''Calculating mean of each type of cell's simulation'''
    TCR_mean[j] = np.mean(TCR, axis=0) / NumCellsTCR
    IN_mean[j] = np.mean(IN, axis=0) / NumCellsIN
    TRN_mean[j] = np.mean(TRN, axis=0) / NumCellsTRN

    # TCRLENGTH = len(TCR_mean)
    # print TCRLENGTH

# f = open('./Results_v0_0916/TCRmean.dat', 'w')
# np.savetxt(f, TCR_mean)
# f.close()

# Taking mean of rows of TCR_mean, IN_mean and TRN_mean. Stored to same variable.
TCR_mean = np.mean(TCR_mean, axis=0) / loop
IN_mean = np.mean(IN_mean, axis=0) / loop
TRN_mean = np.mean(TRN_mean, axis=0) / loop


''' PLOTTING THE MEMBRANE POTENTIAL PLOTS FOR THREE POPULATIONS'''

t = np.arange(0, TotalDataLength, 1)
n_plots = 3
plot = 1

f1 = plt.figure(1)
plt.title('Time series plots for TCR, IN and TNR', fontsize='12')

plt.subplot(n_plots,1,plot)
plot += 1
plt.plot(t, TCR_mean, linewidth='1', linestyle="-", marker='', color='purple', label='TCR')
plt.legend('TCR',fontsize='10')
plt.tick_params(axis='both', which='major', labelsize='10')

plt.subplot(n_plots,1,plot)
plot += 1
plt.plot(t, TRN_mean, linewidth='1', linestyle="-", marker='', color='orange', label='TRN')
plt.legend('TRN',fontsize='10')
plt.tick_params(axis='both', which='major', labelsize='10')

plt.subplot(n_plots,1,plot)
plot += 1
plt.plot(t, IN_mean, linewidth='1', linestyle="-", marker='', color='violet', label='IN')
plt.legend('IN',fontsize='10')
plt.tick_params(axis='both', which='major', labelsize='10')

plt.show()

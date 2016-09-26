''' 
LOG 15th September Basab: Using the code that Akash wrote last year to visualise the data generated with SysManCyberB_Toni_Design_September16.py - looking at the correlation of the spike raster of the three LGN cell populations.

Code Copyright of Akash Bhattacharya, 2nd Year Student, MSc Course in Chemistry with Molecular Physics, Imperial College London, 2015-2016.
'''

#!/usr/bin/python
import numpy as np

import matplotlib.pylab as plt
from pylab import *
from pyNN.random import NumpyRNG, RandomDistribution


def check_if_len_zero(*args):
    verdict = [True] * len(args)
    for i in range(0, len(args)):
        if len(args[i]) > 0:
            verdict[i] = False
    return verdict


'''KNOWN QUANTITIES'''
scale_fact=10
NumCellsTCR = 5*scale_fact ##'''Number of Cells of TCR'''
NumCellsTRN = 4*scale_fact ##'''Number of Cells of TRN'''
NumCellsIN = 1*scale_fact ##'''Number of Cells in IN'''

TimeInt = 10 ##'''millisecond'''
TotalDuration = 2000 ##'''milliseconds; total simulation time'''
TotalDataLength=TotalDuration * TimeInt


''' Creating zeroes-array to store mean values'''
TCR_mean = np.array([0]*TotalDataLength)
TRN_mean = np.array([0]*TotalDataLength)
IN_mean = np.array([0]*TotalDataLength)

'''Load Data'''
# Each simulation of the SysManCyberB_Tonic_Design_September16.py file is run for 5 times. So Loop can be any number from 1 to 5

loop = input('Enter any number from 1 to 20: \n')

#Loading files
spike_source_ex = np.loadtxt('./Sim1_3hz_0916/spikesource_ex_'+`loop`+'.dat')
spike_source_inh = np.loadtxt('./Sim1_3hz_0916/spikesource_inh_'+`loop`+'.dat')
TCRspikes = np.loadtxt('./Sim1_3hz_0916/TCRspikes_'+`loop`+'.dat')
INspikes = np.loadtxt('./Sim1_3hz_0916/INspikes_'+`loop`+'.dat')
TRNspikes = np.loadtxt('./Sim1_3hz_0916/TRNspikes_'+`loop`+'.dat')


TCRmempot = np.loadtxt('./Sim1_3hz_0916/TCRmempot_'+`loop`+'.dat')
INmempot = np.loadtxt('./Sim1_3hz_0916/INmempot_'+`loop`+'.dat')
TRNmempot = np.loadtxt('./Sim1_3hz_0916/TRNmempot_'+`loop`+'.dat')

check_results = check_if_len_zero(spike_source_ex, spike_source_inh, TCRspikes, INspikes, TRNspikes)

'''MATRICES FORMATION'''
signalTCR = []; signalTRN = []; signalIN = []
for i in range(0, len(TCRmempot)):
    signalTCR.append(TCRmempot[i][2])
for i in range(0, len(TRNmempot)):
    signalTRN.append(TRNmempot[i][2])
for i in range(0, len(INmempot)):
    signalIN.append(INmempot[i][2])
signalTCR = np.array(signalTCR)
signalTRN = np.array(signalTRN)
signalIN = np.array(signalIN)
TCR = signalTCR.reshape(NumCellsTCR, TotalDataLength)
IN = signalIN.reshape(NumCellsIN, TotalDataLength)
TRN = signalTRN.reshape(NumCellsTRN, TotalDataLength)

'''Calculating mean of each type of cell's simulation'''

TCR_mean = np.mean(TCR, axis=0)/NumCellsTCR
TCRLENGTH=len(TCR_mean)
print TCRLENGTH

IN_mean = np.mean(IN, axis=0)/NumCellsIN
TRN_mean = np.mean(TRN, axis=0)/NumCellsTRN

#f = open('./Results_v0_0916/TCRmean.dat', 'w')
#np.savetxt(f, TCR_mean)
#f.close()



'''### PLOTTING THE SPIKE SOURCE EXCITATORY AND TCR'''
'''# First column of data is the times the neurons spike.'''
time_sp_ex = spike_source_ex[:, 0]
time_tcrsp = TCRspikes[:, 0]

'''# Second column of data is index of neuron.'''
neuron_sp_ex = spike_source_ex[:, 1]
neuron_tcrsp = TCRspikes[:, 1]


'''# Plotting points, with no line.'''
f1 = plt.figure(1)
plt.subplot(3, 1, 1)
plt.title('Raster plots of TCR spikes (bottom) corresponding to the excitatory source (top)', fontsize='15')
plt.plot(time_sp_ex, neuron_sp_ex, marker='.', linewidth='0', color='orange')
# Scaling y-axis for easy viewing.
plt.ylim(-1.5, neuron_sp_ex[-1] + 1.5)
plt.xlim(-1, TotalDuration)
# Labelling graph.
plt.xlabel("Time/ms")
plt.ylabel("Neuron index")

plt.subplot(3,1,2)
plt.plot(time_tcrsp, neuron_tcrsp, marker='.', linewidth = '0', color='b')
# Scaling y-axis for easy viewing.
plt.ylim(-1.5, neuron_tcrsp[-1] + 1.5)
plt.xlim(-1, TotalDuration)
# Labelling graph.
plt.xlabel("Time/ms")
plt.ylabel("Neuron index")

t = np.arange(0, TotalDataLength, 1)

plt.subplot(3,1,3)
#plt.title('Time series plots for TCR, IN and TNR', fontsize='15')
plt.plot(t, TCR_mean, linewidth='1', linestyle="-", marker='', color='b', label='TCR')
plt.legend('TCR',fontsize='10')
plt.tick_params(axis='both', which='major', labelsize='10')

f1.show()

#'''### PLOTTING THE SPIKE SOURCE INHIBITORY AND IN, TRN'''
## First column of data is the times the neurons spike.
#time_sp_inh = spike_source_inh[:, 0]
#time_insp = INspikes[:, 0]
#if len(TRNspikes) > 0:
#    time_trnsp = TRNspikes[:, 0]
#else:
#    time_trnsp = []
#    print("No spikes from TRN!")

## Second column of data is index of neuron.
#neuron_sp_inh = spike_source_inh[:, 1]
#neuron_insp = INspikes[:, 1]

#if len(TRNspikes) > 0:
#    neuron_trnsp = TRNspikes[:, 1]
#else:
#    neuron_trnsp = []

## Plotting points, with no line.
##f2 = plt.figure(2)
##plt.subplot(3,1,1)
##plt.title('Raster plots of IN and TRN spikes (bottom) corresponding to the inhibitory source (top)', fontsize='15')
##plt.plot(time_sp_inh, neuron_sp_inh, marker = '.', linewidth = '0',color='violet')
### Scaling y-axis for easy viewing.
##plt.ylim(-1.5, neuron_sp_inh[-1] + 1.5)
##plt.xlim(-1, TotalDuration)
### Labelling graph.
##plt.xlabel("Time/ms")
##plt.ylabel("Neuron index")

##plt.subplot(3,1,2)
##plt.plot(time_insp, neuron_insp, marker = '.', linewidth = '0',color='g')
### Scaling y-axis for easy viewing.
##plt.ylim(-1.5, neuron_insp[-1] + 1.5)
##plt.xlim(-1, TotalDuration)
### Labelling graph.
##plt.xlabel("Time/ms")
##plt.ylabel("Neuron index")


##plt.subplot(3,1,3)
##plt.plot(time_trnsp, neuron_trnsp, marker = '.', linewidth = '0',color='cyan')
### Scaling y-axis for easy viewing.
##if len(neuron_trnsp) > 0:
##    plt.ylim(-1.5, neuron_trnsp[-1] + 1.5)
##else:
##    plt.ylim(-0.5, 0.5)

##plt.xlim(-1, TotalDuration)
### Labelling graph.
##plt.xlabel("Time/ms")
##plt.ylabel("Neuron index")

##f2.show()


### VISUALISING THE RECORDED MEMBRANE POTENTIAL
#t = np.arange(0, TotalDataLength, TimeInt)

#f3 = plt.figure(3)
#plt.subplot(3,1,1)
#plt.title('Time series plots for TCR, IN and TNR', fontsize='15')
#plt.plot(t, TCR_mean, linewidth='1', linestyle="-", marker='', label='TCR')
###plt.legend('TCR',fontsize='10')
#plt.tick_params(axis='both', which='major', labelsize='10')

#plt.subplot(3,1,2)
#plt.plot(t, IN_mean, linewidth='1', linestyle="-", color='g', marker='', label='IN')
#plt.ylabel('Membrane Potential (mV)', fontsize='15')
###plt.legend('IN', fontsize='10')
#plt.tick_params(axis='both', which='major', labelsize='10')

#plt.subplot(3,1,3)
#plt.plot(t, TRN_mean, linewidth='1', linestyle="-", marker='', color='cyan', label='TRN')
#plt.xlabel('Time/s', fontsize='15')
###plt.legend('TRN', fontsize='10')
#plt.tick_params(axis='both', which='major', labelsize='10')

#f3.show()
###savefig('./Results/Plot_AVG.png', dpi=1200)

#f4=plt.figure(4)

#t = np.arange(0, TotalDataLength, TimeInt)
#plt.subplot(3,1,1)
#plt.plot(time_trnsp, neuron_trnsp, marker = '.', linewidth = '0',color='cyan', label='TRN')
## Scaling y-axis for easy viewing.
#plt.ylim(-1.5, neuron_trnsp[-1] + 1.5)
#plt.xlim(-1, TotalDuration)
## Labelling graph.
##plt.xlabel("Time/ms")
#plt.ylabel("TRN Neuron index")

#plt.subplot(3,1,2)
#plt.plot(time_tcrsp, neuron_tcrsp, marker = '.', linewidth = '0',color='b', label='TCR')
## Scaling y-axis for easy viewing.
#plt.ylim(-1.5, neuron_tcrsp[-1] + 1.5)
#plt.xlim(-1, TotalDuration)
## Labelling graph.
##plt.xlabel("Time/ms")
#plt.ylabel("TCR Neuron index")

#plt.subplot(3,1,3)
#plt.plot(time_insp, neuron_insp, marker = '.', linewidth = '0',color='g', label='IN')
## Scaling y-axis for easy viewing.
#plt.ylim(-1.5, neuron_insp[-1] + 1.5)
#plt.xlim(-1, TotalDuration)
## Labelling graph.
##plt.xlabel("Time/ms")
#plt.ylabel("IN Neuron index")

#f4.show()

plt.show()

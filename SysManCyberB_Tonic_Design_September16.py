"""
12th - 16th september: 
Trying to record the projections using .getweights - and finally succeeded today with help from Loci - 
Have saved the files with the periodic inputs and revised earlier work. Need to look into comparison with Poisson spike trains now.


9th September: running using spalloc from it 302

19th August: Passing this on to Loci for creaing multiple instances of the main code
modelling the Lateral Geniculate Nucleus (LGN)
LOG: Basab

5th July: running results with possion (aperiodic) spike trains but at same frequencies as with periodic spike trains.
3rd -4th July - running and storing results with regular spike trains at frequencies 23Hz, 19Hz, 15Hz, 11Hz, 8Hz, 4.5Hz, 3.2Hz. The 19Hz frequency
was particularly used at first, to test the several range of parameters, especially varying p1 and p2 (see below for definition) for all 5 combinations
(see at the beginning of the code). Also, a set of Simxa folders contain data when the TRN to TCR feedback connectivity weight was increased to 4.
This was just to test with the set basal value of 2 (as a result of prior experiments on the spiking behaviours). As before, we note that with all other
parameter values as they are, increasing TRN to TCR feedback from 2 to 4 resulted in erratic response from TCR, which intuitively (needs lots of more simulations
to make a statement, and also with a large number of neurons) feels like due to the increased suppression of spiking neurons in the TCR population,
considering the connectivity is also sparse.
30th June - After some initial hiccups, finally could assimilate the data for p1 and p2, and store in Sim 1. But only one instance out of 5 was
'successful' in the sense that the means of the neurons responded till the end.
Now, the values of
29th JUne - running the population network: with both p1 and p2 (probability of connection of TRN And IN respectively, to TCR, to 0.15
11th June: Testing with Single neurons in the circuit.
6th June -  Retested the tonic mode. The offset current values, when there are no synaptic weights, produce fast spking and regular spiking effects with the defined values of abcd parameters.
However, this offset current cannot be negated.
On the other hand, with offset 0, the synaptic weights, there is a conductance based current weighted with the synapses (refer to Andrew's mail) that causes spiking
in the cells, but the effects of fast spiking and regular spiking are not apparent. Unless of course tested with high frequency input spike train. The FS firing can
reproduce a higher frequency at the output that can the RS cells.
1st June 2016: Pape and McCormick 95 - need to determine the ISI.

31st January 2016: The Bursting experiment upto the single mini-column is done - so starting the tonic testing now.
So the first thing to do is change the range of i to (2,3)

24TH JANUARY 2016: FINALISED THE VISUALISATION AT THE END OF THE FILE FOR GENERATING AND SAVING FIGURES 1, 2 AND 3 IN THE PAPER DRAFT
FOR SUBMISSION TO IEEE TRANSACTIONS IN BIOMEDICAL CIRCUITS AND SYSTEMS
THE THREE CONDITIONS ARE:
(1) TONIC BURSTING FOR ALL CELLS WHEN SPIKE SOURCE IS EXCITATORY FOR ALL CELLS, AND WITH MINIMAL ISI TO ELICIT ONE-TO-ONE RESPONSE IN THE CELLS.
THE CONNECTION WEIGHTS ARE SET AS 1.
(2) SAME SPIKE SOURCE INTERVAL, BUT NOW THE PROJECTED SPIKE TO ALL CELLS IS INHIBITORY. THE CONNECTION WEIGHTS ARE INCREASED TO 10.
(3) THE SPIKE SOURCE PROJECTION IS STILL INHIBITORY, CONNECTION WEIGHT IS INCREASED TO 15;
AND THE INTERSPIKE INTERVAL TO ELICIT ONE-TO-ONE FOR GENERATING POST-INHIBITORY-BURST-RESPONSE

LOG JANUARY 2016: SETTING SINGLE NEURON BEHAVIOUR OF TCR, TRN, IN PRIOR TO BUILDING POPULATIONS
LOG 28TH DECEMBER 2015: BEING TESTED WITH THALAMIC FIRING MODES PAPER RESULTS AND OBSERVATIONS: BASAB

LOG: 7TH AUGUST 2015
THE PARAMETETCR OF THE EARLIER MODEL IS CHANGED ACCORDING TO THE ONES WHEN WORKING WITH RAHMI -  MAINLY THE I-OFFSET IS CHANGED TO 0, BECAUSE WHEN THE NEURONS ARE USED
WITHIN A NETWORK, THEIR INPUT CURRENTS WILL BE PROVIDED BY OTHER NEURONS WITHIN THE NETWORK AND HENCE THEY WOULD NOT NEED AN EXPLICIT DRIVER CURRENT.
ALSO, THE U VARIABLES ARE MADE NON-ZERO, WHICH MAKES A DIFFERENCE TO THE INITIAL CONDITION OF THE PLOT BUT AN INSIGNIFICANT CHANGE.
FURTHER, THE TAU PARAMETERS ARE INCREASED, WHICH DOES NOT MAKE ANY VISIBLE DIFFERENCE IN THE CASE OF INHIBITORY SPIKE SOURCE,
BUT IN THE CASE OF AN EXCITATORY SPIKE SOURCE, THE RATE OF SPIKING IN ALL OF THE NEURONS INCREASES.
"""

# !/usr/bin/python

import numpy
import matplotlib.pylab as plt
import pylab
from pylab import *
import pyNN.spiNNaker as p
from pyNN.random import NumpyRNG, RandomDistribution
'''This is the probabilty of connectivity for the TRN'''
p1=[0.25, 0.2, 0.15, 0.11, 0.07]
'''This is the probabilty of connectivity for the IN'''
p2=[0.07, 0.11, 0.15, 0.2, 0.25]



TotalDuration=10000
TimeInt=1
loop=0
while loop<20:
    print ('WE ARE IN LOOP:\n')
    print loop
    
    
    p.setup(timestep=1.0, min_delay=1.0, max_delay=90.0)
    loop += 1
    # Tonic mode parameters
    tcr_a_tonic = 0.02
    tcr_b_tonic = 0.2
    tcr_c_tonic = -65
    tcr_d_tonic = 6
    tcr_v_init_tonic = -65

    in_a_tonic = 0.1
    in_b_tonic = 0.2
    in_c_tonic = -65
    in_d_tonic = 2
    in_v_init_tonic = -70

    trn_a_tonic = 0.02
    trn_b_tonic = 0.2
    trn_c_tonic = -65
    trn_d_tonic = 6
    trn_v_init_tonic = -75

    tcr_a = tcr_a_tonic
    tcr_b = tcr_b_tonic
    tcr_c = tcr_c_tonic
    tcr_d = tcr_d_tonic
    tcr_v_init = tcr_v_init_tonic
    
    
    in_a = in_a_tonic
    in_b = in_b_tonic
    in_c = in_c_tonic
    in_d = in_d_tonic
    in_v_init = in_v_init_tonic
    
    
    trn_a = trn_a_tonic
    trn_b = trn_b_tonic
    trn_c = trn_c_tonic
    trn_d = trn_d_tonic
    trn_v_init = trn_v_init_tonic
   
   
   
    tcr_u_init = tcr_b * tcr_v_init
    in_u_init = in_b * in_v_init
    trn_u_init = trn_b * trn_v_init

    currentPulse = 5
    noCurrentPulse = 0
    # THALAMOCORTICAL RELAY CELLS (TCR)

    TCR_cell_params = {'a': tcr_a_tonic, 'b': tcr_b, 'c': tcr_c, 'd': tcr_d,
                       'v_init': tcr_v_init, 'u_init': tcr_u_init,
                       'tau_syn_E': 10, 'tau_syn_I': 10,
                       'i_offset': noCurrentPulse
                       }

    # THALAMIC INTERNEURONS (IN)

    IN_cell_params = {'a': in_a, 'b': in_b, 'c': in_c, 'd': in_d,
                      'v_init': in_v_init, 'u_init': in_u_init,
                      'tau_syn_E': 10, 'tau_syn_I': 10,
                      'i_offset': noCurrentPulse
                      }

    # THALAMIC RETICULAR NUCLEUS (TRN)

    TRN_cell_params = {'a': trn_a, 'b': trn_b, 'c': trn_c, 'd': trn_d,
                       'v_init': trn_v_init, 'u_init': trn_u_init,
                       'tau_syn_E': 10, 'tau_syn_I': 10,
                       'i_offset': noCurrentPulse
                       }

    '''DEFINING THE POPULATIONS'''
    p.set_number_of_neurons_per_core("IZK_curr_exp", 100)
    scale_fact = 10
    NumCellsTCR = 5*scale_fact
    NumCellsIN = 1*scale_fact
    NumCellsTRN = 4*scale_fact
    TCR_pop = p.Population(NumCellsTCR, p.IZK_curr_exp, TCR_cell_params, label='TCR_pop')
    IN_pop = p.Population(NumCellsIN, p.IZK_curr_exp, IN_cell_params, label='IN_pop')
    TRN_pop = p.Population(NumCellsTRN, p.IZK_curr_exp, TRN_cell_params, label='TRN_pop')



    ''' PERIODIC SPIKE TRAIN INPUT: 23Hz: 44msec isi; 19 Hz: 52 msec isi; 15 Hz:67 msec isi; 11Hz: 91 msec isi; 8 Hz: 125 msec isi; 4.5 Hz: 223 msec isi; 3 Hz: 333msec isi'''

    spike_source_ex = p.Population(NumCellsTCR, p.SpikeSourceArray, {'spike_times': [i for i in range(1500,8000,333)]}, label='spike_source_ex')
    
    
    spike_source_inh = p.Population(NumCellsIN, p.SpikeSourceArray, {'spike_times': [i for i in range(400,5000,333)]}, label='spike_source_inh')


    ''' A-PERIODIC SPIKE TRAIN INPUT'''
#    Rate_Inp = 3
#    spike_source_ex = p.Population(NumCellsTCR, p.SpikeSourcePoisson, {'rate':Rate_Inp, 'duration':840,'start':150}, label='spike_source_ex')
#    spike_source_inh = p.Population(NumCellsIN, p.SpikeSourcePoisson, {'rate': Rate_Inp, 'duration':250,'start':250}, label='spike_source_inh')


    # DEFINING THE PROJECTION AND WEIGHT PARAMETERS
    projList=list()
    
    tcr_weights = 6
    in_weights = 2
    
    
    
    '''Source2TCR'''
    Proj0 = p.Projection(spike_source_ex, TCR_pop, p.FixedProbabilityConnector(p_connect=0.07, weights=tcr_weights, delays=5), target='excitatory')
    projList.append(Proj0)
    
    '''Source2IN'''
    Proj1 = p.Projection(spike_source_inh, IN_pop, p.FixedProbabilityConnector(p_connect=0.47, weights=in_weights, delays=5), target='excitatory')
    projList.append(Proj1)


    '''TCR2TRN'''
    Proj2 = p.Projection(TCR_pop, TRN_pop, p.FixedProbabilityConnector(p_connect=0.35, weights=4, delays=3), target='excitatory')
    projList.append(Proj2)
    
    '''TRN2TCR''' 
    Proj3 = p.Projection(TRN_pop, TCR_pop, p.FixedProbabilityConnector(p_connect=0.15, weights=2, delays=3), target='inhibitory')
    projList.append(Proj3)
    
    '''TRN2TRN'''
    Proj4 = p.Projection(TRN_pop, TRN_pop, p.FixedProbabilityConnector(p_connect=0.15, weights=2, delays=1), target='inhibitory')
    projList.append(Proj4)


    '''IN2TCR'''
    Proj5 = p.Projection(IN_pop, TCR_pop, p.FixedProbabilityConnector(p_connect=0.16, weights=4, delays=2), target='inhibitory')
    projList.append(Proj5)
    
    '''IN2IN'''
    Proj6  = p.Projection(IN_pop, IN_pop, p.FixedProbabilityConnector(p_connect=0.24, weights=2, delays=1), target='inhibitory')
    projList.append(Proj6)

    # RECORDING THE MEMBRANE POTENTIAL
    
    ''' RECORDING THE MEMBRANE POTENTIAL'''

    TCR_pop.record_v()
    IN_pop.record_v()
    TRN_pop.record_v()

    ''' RECORDING THE SPIKES'''
    spike_source_ex.record()
    spike_source_inh.record()
	

    TCR_pop.record()
    IN_pop.record()
    TRN_pop.record()


    p.run(TotalDuration)

    ''' STORING DATA IN TEXT FILES'''



    print ('SAVING THE FILES')
    TCR_pop.print_v('./Sim3_3hz_0916/TCRmempot_'+`loop`+'.dat')
    IN_pop.print_v('./Sim3_3hz_0916/INmempot_'+`loop`+'.dat')
    TRN_pop.print_v('./Sim3_3hz_0916/TRNmempot_'+`loop`+'.dat')
    spike_source_ex.printSpikes('./Sim3_3hz_0916/spikesource_ex_'+`loop`+'.dat')
    spike_source_inh.printSpikes('./Sim3_3hz_0916/spikesource_inh_'+`loop`+'.dat')
    TCR_pop.printSpikes('./Sim3_3hz_0916/TCRspikes_'+`loop`+'.dat')
    IN_pop.printSpikes('./Sim3_3hz_0916/INspikes_'+`loop`+'.dat')
    TRN_pop.printSpikes('./Sim3_3hz_0916/TRNspikes_'+`loop`+'.dat')
    print('SAVED THE FILES \n \n \n')
    

    ''' 
    PROJECTION WEIGHTS RECORDED IN A .CSV FILE WHICH CAN BE OPENED IN EXCEL.
    FOR PLOTTING IN MATLAB, JUST TYPE LOAD <FILENAME>.CSV WHEN THIS SHOULD SHOW UP IN WORKSPACE.
    USING THE PLOT TOOLBAR IN LATER VERSIONS OF MATLAB, JUST CLICKING ON 'SPY' SHOULD GIVE A SCATTER 
    PLOT OF THE 2D MATRIX - EVEN THOUGH THE NON-CONNECTED ELEMENTS ARE AS NaN (NOT A NUMBER). 
    HOWEVER, A SIMPLE COMMAND LIKE FIND() IN MATLAB CAN IDENTIFY THE NON-NUMERIC ELEMENTS AND SET THEM 
    TO ZERO.
    '''
    counter=0;
    for projLoop in projList:
    	nowArr = projLoop.getWeights(format='array', gather=True)
    	counter=counter+1
    	np.savetxt('thisloop'+`loop`+'_thisproj'+`counter`+'.csv', nowArr, delimiter=',')
    	
    	
    
    p.end()




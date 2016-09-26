"""

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
while loop<1:
    print ('WE ARE IN LOOP:\n')
    print loop
    
    
    p.setup(timestep=0.1, min_delay=1.0, max_delay=14.0)
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



    ''' PERIODIC SPIKE TRAIN INPUT: 23Hz: 44msec isi; 19.2 Hz: 52 msec isi; 15 Hz:67 msec isi; 11Hz: 91 msec isi; 8 Hz: 125 msec isi; 4.5 Hz: 223 msec isi; 3 Hz: 333msec isi'''

    spike_source_ex = p.Population(NumCellsTCR, p.SpikeSourceArray, {'spike_times': [i for i in range(1500,8000,125)]}, label='spike_source_ex')
    
    
    spike_source_inh = p.Population(NumCellsIN, p.SpikeSourceArray, {'spike_times': [i for i in range(500,5000,125)]}, label='spike_source_inh')


    ''' A-PERIODIC SPIKE TRAIN INPUT'''
#    Rate_Inp = 3
#    spike_source_ex = p.Population(NumCellsTCR, p.SpikeSourcePoisson, {'rate':Rate_Inp, 'duration':6500,'start':1500}, label='spike_source_ex')
#    spike_source_inh = p.Population(NumCellsIN, p.SpikeSourcePoisson, {'rate': Rate_Inp, 'duration':4600,'start':400}, label='spike_source_inh')


    # DEFINING THE PROJECTION AND WEIGHT PARAMETERS
    projList=list()
    
    tcr_weights = 6
    in_weights = 6
    
    
    
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
    Proj3 = p.Projection(TRN_pop, TCR_pop, p.FixedProbabilityConnector(p_connect=0.07, weights=2, delays=3), target='inhibitory')
    projList.append(Proj3)
    
    '''TRN2TRN'''
    Proj4 = p.Projection(TRN_pop, TRN_pop, p.FixedProbabilityConnector(p_connect=0.15, weights=2, delays=1), target='inhibitory')
    projList.append(Proj4)


    '''IN2TCR'''
    Proj5 = p.Projection(IN_pop, TCR_pop, p.FixedProbabilityConnector(p_connect=0.24, weights=6, delays=2), target='inhibitory')
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
    n= 1  #input('enter number of simulation')
    f= 8   #input('enter frequency of simulation')
    foldername="Sim%d_%dhz_0916" % (n,f)
    TCR_pop.print_v('./'+foldername+'/TCRmempot_'+`loop`+'.dat')
    IN_pop.print_v('./'+foldername+'/INmempot_'+`loop`+'.dat')
    TRN_pop.print_v('./'+foldername+'/TRNmempot_'+`loop`+'.dat')
    spike_source_ex.printSpikes('./'+foldername+'/spikesource_ex_'+`loop`+'.dat')
    spike_source_inh.printSpikes('./'+foldername+'/spikesource_inh_'+`loop`+'.dat')
    TCR_pop.printSpikes('./'+foldername+'/TCRspikes_'+`loop`+'.dat')
    IN_pop.printSpikes('./'+foldername+'/INspikes_'+`loop`+'.dat')
    TRN_pop.printSpikes('./'+foldername+'/TRNspikes_'+`loop`+'.dat')
    print('SAVED THE FILES \n \n \n')
    

    ''' 
    PROJECTION WEIGHTS RECORDED IN A .CSV FILE WHICH CAN BE OPENED IN EXCEL.
    FOR PLOTTING IN MATLAB, JUST TYPE LOAD <FILENAME>.CSV WHEN THIS SHOULD SHOW UP I#    print ('SAVING THE FILES')
#    n=input('enter number of simulation')
#    f=input('enter frequency of simulation')
#    foldername="Sim%d_%dhz" % (n,f)
#    TCR_pop.print_v('./'+foldername+'/TCRmempot_'+`loop`+'.dat')N WORKSPACE.
    USING THE PLOT TOOLBAR IN LATER VERSIONS OF MATLAB, JUST CLICKING ON 'SPY' SHOULD GIVE A SCATTER 
    PLOT OF THE 2D MATRIX - EVEN THOUGH THE NON-CONNECTED ELEMENTS ARE AS NaN (NOT A NUMBER). 
    HOWEVER, A SIMPLE COMMAND LIKE FIND() IN MATLAB CAN IDENTIFY THE NON-NUMERIC ELEMENTS AND SET THEM 
    TO ZERO.
    '''
#    counter=0;
#    for projLoop in projList:
#    	nowArr = projLoop.getWeights(format='array', gather=True)
#    	counter=counter+1
#    	np.savetxt('thisloop'+`loop`+'_thisproj'+`counter`+'.csv', nowArr, delimiter=',')
    	
    	
    
    p.end()




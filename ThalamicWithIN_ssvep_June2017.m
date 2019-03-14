%% LOG JUNE 2017, BASAB: USING THIS TO GENERATE RESULTS CORRESPONDING TO PERIODIC+
% NOISY INPUT - WITH THE PURPOSE TO GENERATE HERMANN STYLE FIGURE WITH
% OUTPUT FROM THIS MODEL - LIKE WE DID FOR THE SPINNAKER BASED MODEL. SO WE
% JUST SAVE THE POWER SPECTRA OF THE TCR FROM HERE AND THEN WILL CALL THIS
% FROM WITHIN THE HERMANN STYLE PLOT SCRIPT.

% WITH IMPULSE INPUT STRENGTH 10, THE POWER SPECTRA SHOWS A STRONG
% EXISTENCE OF LOWER SPECTRA COMPONENTS FOR ALL FREQUENCIES, AND NOT MUCH
% HARMONIC CONTENTS.
% WITH IMPULSE INPUT STRENGTH 25, THE POWER SPECTRA SHOWS HARMONICS.

% NEXT WE REMOVE THE IN.
% AND FIRST TEST WITH IMPULSE INPUT STRENGTH 25.

%% LOG OCTOBER 2016: 4TH REVISION OF PAPER - AND REVISITING THE RESULTS 
% CORRESPONDING TO VARYING CONNECTIVITY PARAMETERS - TABLE 3 IN PAPER. MADE
% A FEW CHANGES TO VISUALISATION AND GENERATION OF DATA: (1) SIMULATION
% TIME IS 40 SECONDS WITH CUT OFF AT 10 AND 30 SECONDS.
%% LOG MARCH 2016: FIRST REVISION OF THE FRONTIERS IN COMPUTATIONAL
%% NEUROSCIENCE PAPER
% THE CODE IS COPYRIGHT OF BASABDATTA SEN BHATTACHARYA; PLEASE ACKNOWLEDGE THE
% AUTHOR IF USING THE CODE IN PART OR FULL.


tic
clear all
close all

%% MEMBRANE CAPACITANCE IN MICROFARADS
Cm=1;

%% THE TIME PARAMETERS AND VECTOR
delt=0.001; %% 1 millisecond
endtime=40; %% seconds

timevec=0:delt:endtime;
timelen=numel(timevec);

%% THE STEPS OF SOLUTION
mu=delt;

%% GENERATING NOISY RETINAL INPUT
V_ret_mean=-65;
V_ret_std=2;

%% TRANSMITTER CONCENTRATION PARAMETERS INITIALISED
Tmax=1; %% in Moles (M): so 2.84mM; from fig. 5 legend 1994 destexhe
Kp = 3.8; %% in mV
Vp=-32; %% in mV
Tval=[Tmax Kp Vp];
%% AMPA PARAMETERS INITIALISED
E_ampa=0;

alpha_ampa=1000; %% mM^(-1)ms^(-1)
beta_ampa=50;%% ms^(-1)%%
g_ampa_ret2tcr=300; 
g_ampa_ret2in=100; 
g_ampa_tcr2trn=100; 
ampaval=[E_ampa alpha_ampa beta_ampa g_ampa_tcr2trn g_ampa_ret2tcr g_ampa_ret2in];

%% GABA_A PARAMETERS INITIALISED
alpha_gaba_a=1000; %% mM^(-1)ms^(-1)
beta_gaba_a=40;%% ms^(-1)
g_gaba_a_trn2tcr=100;
g_gaba_a_trn2trn=100;
E_gaba_a_trn2tcr=-85; %%
E_gaba_a_trn2trn=-75;
g_gaba_a_in2tcr=g_gaba_a_trn2tcr;
g_gaba_a_in2in=100;
E_gaba_a_in2tcr=E_gaba_a_trn2tcr; %%
E_gaba_a_in2in=E_gaba_a_trn2trn;
gabaAval=[alpha_gaba_a beta_gaba_a g_gaba_a_trn2tcr g_gaba_a_trn2trn E_gaba_a_trn2tcr E_gaba_a_trn2trn g_gaba_a_in2tcr g_gaba_a_in2in E_gaba_a_in2tcr E_gaba_a_in2in];

%% GABA_B PARAMETERS INITIALISED
alpha1_gaba_b=10; %% mM^(-1)ms^(-1)
beta1_gaba_b=25;%% ms^(-1)
alpha2_gaba_b=15; %% mM^(-1)ms^(-1)
beta2_gaba_b=5;%% ms^(-1)
g_gaba_b =60;%% in milliSiemens (mS)

E_gaba_b=-100;
Kd_gaba_b=100;
n=4;

gabaBval=[alpha1_gaba_b  beta1_gaba_b  alpha2_gaba_b  beta2_gaba_b   g_gaba_b   E_gaba_b   Kd_gaba_b   n];
%% LEAK PARAMETERS INITIALISED
E_leak_tcr=-55; % Leak voltage
g_leak_tcr=10; %%

E_leak_trn=-72.5; % Leak voltage
g_leak_trn=10;
% g_leak_trn_arr=[0.0025 0.0075 0.025 0.075 0.25 0.75];

E_leak_in=-72.5; % Leak voltage
g_leak_in=10;

leakval=[E_leak_tcr g_leak_tcr E_leak_trn g_leak_trn E_leak_in g_leak_in];
%% VECTOR INITIALISATION VALUES
% Ginit=1;
Ginit=0.001;

r_ampa_initval=Ginit;
r_gaba_a_initval=Ginit;
r_gaba_b_initval=Ginit;

vinit_tcr=-65; %% in mV
vinit_trn=-85; %% in mV
vinit_in=-75; %% in mV

initval=[Ginit vinit_tcr vinit_trn vinit_in];

% %% CONNECTIVITY PARAMETERS INITIALISED

Cnte=35;
Ctnia=(3/8)*30.9;%% Half of the the inhibitory synaptic connections from the RE pop to the TCR pop
Ctnib=(1/8)*30.9;
Ctre=7.1;
Cnsi=20;
Cire = 47.4;
Cisi = 23.6;
%% THIS IS THE PARAMETER THAT WE WILL BE DELETING
Ctii = 15.45;

Carr=[Cnte Ctnia Ctnib Ctre Cnsi Cire Ctii Cisi];

%% EVOLUTION OF THE MEMBRANE POTENTIAL: TCR and TRN
freqrange = 10:1:12;%1:1:50;
Vtcravgmat = zeros(length(freqrange),timelen);
for freq = 1:length(freqrange)
    int_hz=round(1000/freqrange(freq));
    noftrials=20;
    startind=9001; endind=39001;
    
    Vtcrmat=zeros(noftrials,timelen);
    
    
    %% FREQUENCY ANALYSIS
    Fs = 1000;
    NFFT=4*Fs;
    WindowType = 'hamming';
    SegmentLength=(1/4)*Fs;
    OverlapPercent=50;
    Normalised=0;
    hp = spectrum.welch(WindowType,SegmentLength,OverlapPercent);
    
    
    %% Construct an FDESIGN object and call its BUTTER method.
    Fc1 = 1; % First Cutoff Frequency
    Fc2 = 200; % Second Cutoff Frequency
    N = 10; % Order
    h = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
    Hd = design(h, 'butter');
    [B,A]=sos2tf(Hd.sosMatrix,Hd.Scalevalues);
    
    fprintf('At the %d-th frequency \n',freqrange(freq))
    for numtrial=1:noftrials
        fprintf('within the %d-th loop \n',numtrial)
        
        %% TRANSMITTER CONCENTRATION VECTOR CORRESPONDING TO RETINAL NOISY INPUT
        R=randn(1,timelen);
        V_ret = ((((R-mean(R)) ./ std(R)) .* V_ret_std) + V_ret_mean);%V_ret_mean.*ones(1,timelen);%
        %     T_ret=Tmax./(1+exp(-(V_ret-Vp)./Kp));
        
        %% EVENT INPUTS AT DIFFERENT FREQUENCIES
        Vimpulse = zeros(1, length(V_ret));
        % ************************************************************************
        Vimpulse(int_hz:int_hz:end)=20;%%   STRENGHT OF INPUT %% STRENGTH OF INPUT <<<<<<<<<<<<<---*****************
        % ***************************************************************************
        V_eventinp = V_ret + Vimpulse;
        
        
        T_ret=Tmax./(1+exp(-(V_eventinp-Vp)./Kp));
        
        
        Y=rk45func_thalmodwithin(Carr, initval, Tval, ampaval, gabaAval, gabaBval, leakval, T_ret, Cm);
        
        Vtcrmat(numtrial,:)=Y(10,1:end-1);
                  
        
        %% Filter data
        filtData1 = filtfilt(B,A,Vtcrmat(numtrial,startind:endind));
        
        hpopts1 = psdopts(hp,filtData1);
                
        set(hpopts1,'Fs',Fs,'NFFT',NFFT,'Normalized',Normalised)
                
        hpsd1 = psd(hp,filtData1,hpopts1);
                
        Ptcrmat(numtrial,:)=hpsd1.Data';
        
    end
    Vtcravgmat(freq,:)=mean(Vtcrmat,1);
    
    Ptcravgmat(freq,:)=mean(Ptcrmat,1);
    
end
fr=hpsd1.Frequencies;
% save thalNoin_ssvep_1_50_inp10.mat Vtcravgmat Ptcravgmat
toc


figure,imagesc(1:50,fr(4:201),Ptcravgmat(:,4:201)'),axis('xy'), colormap('hsv'),colorbar

figure,plot(timevec,V_eventinp),box on, title('MODEL INPUT')
figure, plot(fr,Ptcravgmat,'r'),box on, title('MODEL OUTPUT POWER SPECTRA')

figure,plot(timevec(15000:20000), Vtcravgmat(1,15000:20000),'c')
hold on, plot(timevec(15000:20000),Vtcravgmat(2,15000:20000),'g')
hold on, plot(timevec(15000:20000),Vtcravgmat(3,15000:20000),'m')
title('OUTPUT OF THE MODEL')


%% Visualising power spectra
% THE TIME SERIES
% figure,subplot(3,1,1), plot(timevec(startind:endind),Vtcravgmat(1,startind:endind),'r')%ylabel('TCR','fontsize',14)
% xlabel('Time (seconds)','fontsize',14)
% 
% set(gca,'Fontsize',12),box off
% 
% hold on, subplot(3,1,2), plot(timevec(startind:endind),Vinavgmat(1,startind:endind),'m')%ylabel('TCR','fontsize',14)
% xlabel('Time (seconds)','fontsize',14)
% 
% set(gca,'Fontsize',12), box off
% 
% hold on, subplot(3,1,3), plot(timevec(startind:endind),Vtrnavgmat(1,startind:endind),'b') %ylabel('Membrane Potential (mV)','fontsize',14)
% xlabel('Time (seconds)','fontsize',14)
% 
% set(gca,'Fontsize',12),box off

% THE FOURIER TRANSFORM AND POWER SPECTRAL DENSITY
% fr=hpsd1.Frequencies;
% 
%     figure, subplot (3,1,1), plot(fr,Ptcravg,'r'),box on, title('TCR')
%     xlim([0 50])
%     
%     set(gca,'Fontsize',12),
%     
%     hold on, subplot(3,1,2), plot(fr,Pinavg,'m'), box on, title('IN')
% %     xlabel('Frequency (Hz)','fontsize',14), 
% %     ylabel('Power spectral density','fontsize',14),
%     xlim([0 50])
%     
%     set(gca,'Fontsize',12),
%     
%     hold on, subplot(3,1,3), plot(fr,Ptrnavg,'b'),box on, title('TRN')
%     xlim([0 50])
%     axis('xy')
%     set(gca,'Fontsize',12),


 
   %% SAVE THE GENERATED DATA AS A .MAT FILE (CASE: WITH IN)  
%    relpowarr_basal=[relpowarr_tcr';relpowarr_in';relpowarr_trn'];
%    Pavgmat_basal=[Ptcravg;Pinavg;Ptrnavg];
%    save basal_output.mat Pavgmat_basal relpowarr_basal
   
   

%% LOG: BASAB 16TH SEPTEMBER BASAB:
% TESTING WITH A NEW SET OF SIMULATION BEING RUN ON THE SPINNAKER SERVER

%LOG 3rd July: The SimX (X=1 to 5) folders are with varying values of P for
% IN to TCR and TRN to TCR connectivity i.e. the inhibitory efferents to
% the TCR population. However, the connectivity strength is set to 4 for IN to
% TCR and 2 for TRN to TCR (i.e. less than the feed-forward pathway from
% the IN). The SimXa folders are for comparing with the case when the
% feed-back pathway connectivity (From the TRN) strenth is also increased
% to 4 from 2; i.e. both inhibitory afferents of the TCR have a
% connectivity strength of 4.
%% LOG: BASAB 30TH JUNE 2016
%% THIS IS THE MAIN MATLAB CODE FOR VISUALISING THE DATA COLLECTED FROM
%%RUNNING SIMULATIONS ON SPINNAKER FOR THE SYSTEMS CYBER MAN B PAPER.
%%(FOR SIMULATION ON SPINNAKER, FIRST THE
%%SYSTEMSMANCYBERMANB_TONIC_DESIGN_JUNE16.PY IS EXECUTED. THE DATA
%%COLLECTED IN THE RESULTS_01062016 FOLDER ARE INVOKED FROM THE
%%REMOVEHEADERLINES.PY FILE AND THE HEADER LINES ARE REMOVED, TO STORE IN
%%THE SIM1 FOLDER WITHIN THE SYSTEMSCYBERMANB FOLDER. ALL OF THESE FILES
%%AND FOLDERS ARE WITHIN THE JUNE2016 FOLDER, WHERE THE CURRENT FILE IS.)
%% THE POPULATIONS ARE FIRST DEFINED AND THE DATA LOADED FROM THE .DAT FILE
%%IS STRUCTURED.


%% Number of neurons in each cell population:
clear all
clc
% close all

scale_fact=10;

% tcrpop=5*scale_fact;
% inpop=1*scale_fact;
trnpop=4*scale_fact;


TotalDuration=10000; % Total duration of simulation is 1000 msec
TimeInt=1/0.1; %%Samplint time is 0.1 millisecond
TotalTimeSteps = TotalDuration * TimeInt; %Total time steps is 10000


%Setting the low and high cut off points. These are divided into two parts
%based on the two intervals when the retinal inputs are applied to the IN
%cell population. The exact times for inhibitory inputs are: 2.5 - 5 sec;
%and 7.5 - 8.5 seconds. The TCR cells receive the retinal inputs from 0.5
%sec to 9.9 sec


locut=100; 
hicut=99900; 

% current_neuronpop=tcrpop;
% current_neuronpop=inpop;    
current_neuronpop=trnpop;



display ('Simulation 8 Hz')
f=8; n=1;
for loop=1:1
   
    
%     filename=sprintf('./SysManCyberB_Sep16/PeriodicSpikeTrain/exp3hz/Sim1/TCRmempot_%d.dat',loop);
%         filename=sprintf('./SysManCyberB_Sep16/PeriodicSpikeTrain/exp3hz/Sim1/INmempot_%d.dat',loop);
%         filename=sprintf('./SysManCyberB_Sep16/PeriodicSpikeTrain/exp3hz/Sim1/TRNmempot_%d.dat',loop);
%   filename=sprintf('../Sim%d_%dhz_0916/TCRmempot_%d.dat',n,f,loop); 
% filename=sprintf('../Sim%d_%dhz_0916/INmempot_%d.dat',n,f,loop);
filename=sprintf('../Sim%d_%dhz_0916/TRNmempot_%d.dat',n,f,loop);
    
    fid = fopen(filename);
    neuron_par = textscan(fid, '%f %f %f'); %% for the rest of the files,
    fclose(fid);
    
    startind=1;
    mempot=zeros(current_neuronpop,TotalTimeSteps);
    for i = 1:current_neuronpop
        try
            mempot(i,:) = neuron_par{1,3}(startind:TotalTimeSteps+startind-1);
            %         mempot(i,:) = neuron_par{1,1}(startind:TotalTimeSteps+startind-1);
            startind=TotalTimeSteps+startind;
        catch
            display('fault at i='),i,display('@'),loop
            
            break
        end
    end
    
    
    %% Visualise Data
    M(loop,:)=mean(mempot,1);
%     M=mean(mempot,1);
    
% end
% M(find(M<-100))=-100;

% figure,plot(100:TotalTimeSteps-100,M(:,100:TotalTimeSteps-100))101
%% visualizing the behaviour of each neuron in the population
% locut=[200 2250 4150 6150 8150];
% hicut=[2200 4100 6100 8100 9950];
% mempot(find(mempot<-100))=-100;
% for len=1:length(locut)
%     figure,imagesc(locut(len):hicut(len),1:current_neuronpop,mempot1(1:current_neuronpop,locut(len):hicut(len)))
end
meanM=mean(M,1);
figure, plot(locut:TimeInt:hicut,M(:,locut:TimeInt:hicut))
hold on, plot(locut:TimeInt:hicut,meanM(1,locut:TimeInt:hicut),'--r','linewidth',0.5)

Fs = 1000;
N   = 10;  % Order
NFFT=4*Fs;
WindowType = 'hamming';
SegmentLength=(1/4)*Fs;
OverlapPercent=50;
Normalised=0;
hp = spectrum.welch(WindowType,SegmentLength,OverlapPercent);

% Construct an FDESIGN object and call its BUTTER method.
Fc1 = 1;   % First Cutoff Frequency
Fc2 = 100;  % Second Cutoff Frequency
h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
Hd = design(h, 'butter');
[B,A]=sos2tf(Hd.sosMatrix,Hd.Scalevalues);

% % figure,
% %% THIS BIT OF THE CODE IS TO SELECTIVELY TEST THE DATA COLLECTED FROM
% %%SPINNAKER FOR ITS POWER SPECTRA BEHAVIOUR; ALL NEURONS DO NOT REACT TO
% %%THE INPUTS AND JUST FALLS INTO A LOW VOLTAGE NOISY MODE. LOG: BASAB 30TH
% %%JUNE 2016
% 
% for ind = 1:size(M,1) % this is the number of rows in M, which is the number of trials on SpiNNaker, and called udn
% %     for len=1:length(locut)
%         filtData(ind,:) = filtfilt(B,A,M(ind,locut:hicut));
% %         filtData(filtData<-3)=-3;
% %         filtData(filtData>5)=5;
%         figure,plot(locut:hicut,filtData(ind,:))
%         pause
%         
%         hpopts = psdopts(hp,filtData(ind,:));
%         set(hpopts,'Fs',Fs,'NFFT',NFFT,'Normalized',Normalised)
%         hpsd = psd(hp,filtData(ind,:),hpopts);
%         Pmat=hpsd.Data';
%         fr=hpsd.Frequencies;
%         figure, plot(fr(1:800),Pmat(1:800),'b','linewidth',2)
%         xlim([-0.5 201])
%         set(gca,'XTick',[1 20 40 80 100],...
%             'XTickLabel','|20|40|80|100', 'Fontsize',12)
%         xlabel('Frequency (Hz)','Fontsize',14)
%         ylabel('Power spectra magnitude','Fontsize',14)
%         pause
%         x=Fc1; y=Fc2;
%         [fr vismat]=fun_stft(M(ind,locut:hicut), x, y);
%         figure, imagesc([],fr((4*x+1):(4*y+1)),vismat((4*x+1):(4*y+1),2:end));
%         xlabel('Time windows','Fontsize',14);
%         ylabel('frequency(Hz)','Fontsize',14);
%         axis('xy')
%         set(gca,'Fontsize',12),
%         ylim([1 200])
%         title('input frequency is variable Hz')
%         pause
% %     end
% %     close all
%    
% end
% 
% % % %************************************************
%% THIS IS THE ORIGINAL CODE TO LOOK INTO THE POWER SPECTRA AS A WHOLE AND
%%ALSO FOR THE WHOLE MATRIX.
    for ms=1:size(M,1)
        filtData(ms,:) = filtfilt(B,A,M(ms,locut:TimeInt:hicut));
        hpopts = psdopts(hp,filtData(ms,:));
        set(hpopts,'Fs',Fs,'NFFT',NFFT,'Normalized',Normalised)
        hpsd = psd(hp,filtData(ms,:),hpopts);
        Pmat(ms,:)=hpsd.Data';
    end
    meanP=mean(Pmat,1);
    fr=hpsd.Frequencies;
    figure, 
    plot(fr(1:400),Pmat(:,1:400))
    hold on, 
    plot(fr(1:400),meanP(1:400),'--g','linewidth',2)
    xlim([-0.5 101])
%     set(gca,'XTick',[1 20 40 80 100],...
%         'XTickLabel','|20|40|80|100', 'Fontsize',12)
    %     'XTickLabel','1||10|20|25|30|40|50|75|100|200', 'Fontsize',12)

    xlabel('Frequency (Hz)','Fontsize',14)
    ylabel('Power spectra magnitude','Fontsize',14)

    mean_filtdata = mean(filtData,1);
    figure, plot (mean_filtdata)
%     plot(locut:hicut,filtData(:,1:end))
    
% %     hold on, 
% %     plot(locut:hicut,mean_filtdata,'g','linewidth',1)

    % STFT
        x=Fc1; y=Fc2;
        [fr, vismat]=fun_stft(meanM(1,locut:TimeInt:hicut), x, y);
        figure, imagesc([],fr((4*x+1):(4*y+1)),vismat((4*x+1):(4*y+1),2:end));
        xlabel('Time windows','Fontsize',14);
        ylabel('frequency(Hz)','Fontsize',14);
        axis('xy')
        set(gca,'Fontsize',12),
        ylim([x y])
        title('input frequency is variable Hz')

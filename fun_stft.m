function [fr, vismat]=fun_stft(Vm, Fc1, Fc2)
 
        fs=1000;
        N   = 10;  % Order
        h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, fs);
        Hd = design(h, 'butter');
        
        
        [B,A]=sos2tf(Hd.sosMatrix,Hd.Scalevalues);
        Data = filtfilt(B,A,Vm);
              
        %% STFT
        
        time_window_len=1000;%% 
        slider_span=0.50*time_window_len; %% window is slided by 125___ elements which is every 125___ millisecond 
        x=size(Data);
%         no_of_windows = round(x(2)/slider_span);
        nfft=4*fs;
        xind = ceil(x(2)/slider_span);
        yind = (nfft/2)+1;
        stft_mat = zeros(xind,yind);
        zeropadding_len=time_window_len - slider_span;
        X = [zeros(1,zeropadding_len) Data zeros(1,zeropadding_len)]; % padding the signal with zeros
        % Hamming window eases out the ripples compared to a rectangular window
        WIN = transpose(hamming(time_window_len)); %ones(1,time_window_len);
        
        iter = 0;
        for i = 1:slider_span:x(2)-1
            iter = iter + 1;
            time_window = X(i:(i + time_window_len - 1));  % making window
            
            signal_window = time_window.* WIN; %   y(t)=h(t)*w(n-t)
            signal_out(iter,:) = abs(fft(signal_window, nfft)); %fft of y(t) gives us stft
            stft_mat(iter,:) = signal_out(iter,1:yind);
        end
       fr = 0:0.25:fs/2;
       vismat=stft_mat';
        
        
function [berVal,throughput] = fn_OFDM_Basic2(subcrs,cyc,sigtonoise,delaySp,modOrder)

NgType=1; % NgType=1/2 for cyclic prefix/zero padding
Ch=1;  % Ch=0/1 for AWGN/multipath channel

if NgType==1
    nt='CP';
elseif NgType==2
    nt='ZP';
end

if Ch==0
    chType='AWGN';
    Target_neb=500;
else
    chType='CH';
    Target_neb=500;
end

% PowerdB=[0 -8 -17 -21 -25]; % Channel tap power profile 'dB'
PowerdB = delaySp;
Delay = [0 3 5 7 8];        % Channel delay 'sample'                        % Discrete path delay in seconds.
Power=10.^(PowerdB/10);     % Channel tap power profile 'linear scale'
Ntap=length(PowerdB);       % Chanel tap number
Lch=Delay(end)+1;           % Channel length

% rayleighChan = comm.RayleighChannel(...
%     'SampleRate',10e3, ...
%     'PathDelays',[0 1.5].*1e-4, ...
%     'AveragePathGains',[1 2],...
%     'NormalizePathGains',true,...
%     'MaximumDopplerShift',30,...
%     'DopplerSpectrum',{doppler('Gaussian',0.6),doppler('Flat')},...
%     'RandomStream','mt19937ar with seed',...
%     'Seed',22,...
%     'PathGainsOutputPort',true);
    

Nbps=modOrder; 
M=2^Nbps;                   % Modulation order=2/4/6 for QPSK/16QAM/64QAM
Nfft=subcrs;                % FFT size                    
Ng=cyc;                     % Guard interval length
Nsym=Nfft+Ng;               % Symbol duration
Nvc=Nfft/4;                 % Nvc=0: no virtual carrier
Nused=Nfft-Nvc;

EbN0=sigtonoise;           % EbN0
N_iter=1e3;       % Number of iterations for each EbN0
Nframe=3;         % Number of symbols per frame
sigPow=0;         % Signal power initialization

throughput = (((Nfft/Nsym)* Nused * Nbps * (3/4))/(4e-6))/1e6; 

% file_name=['OFDM_BER_' chType '_' nt '_' 'GL' num2str(Ng) '.dat'];
% fid=fopen(file_name, 'w+');

norms=[1 sqrt(2) 0 sqrt(10) 0 sqrt(42)];     % BPSK 4-QAM 16-QAM

for i=0:length(EbN0)
   randn('state',0);
   rand('state',0); 
   Ber = 0;
%    Ber2=Ber(); % BER initialization  
   
   Neb=0;
   Ntb=0; % Initialize the number of error/total bits
   
   for m=1:N_iter
      % Tx______________________________________________________________
      X= randi(1,Nused*Nframe,M); % bit: integer vector
      Xmod= qammod(X,M,'gray')/norms(Nbps);
      
      if NgType~=2              % BURASI SIKINTILI GİBİ
          x_GI = zeros(1,Nframe*Nsym);
       elseif NgType==2
          x_GI = zeros(1,Nframe*Nsym+Ng);
        % Extend an OFDM symbol by Ng zeros 
      end
      
      kk1=[1:Nused/2]; 
      kk2=[Nused/2+1:Nused];
      kk3=1:Nfft; 
      kk4=1:Nsym;
      
      for k=1:Nframe
         if Nvc~=0
             X_shift= [0 Xmod(kk2) zeros(1,Nvc-1) Xmod(kk1)];
         else
             X_shift= [Xmod(kk2) Xmod(kk1)];
         end
         
         x= ifft(X_shift);
         x_GI(kk4)= guard_interval(Ng,Nfft,NgType,x);
         kk1=kk1+Nused; 
         kk2=kk2+Nused;
         kk3=kk3+Nfft; 
         kk4=kk4+Nsym;
      end
      
      if Ch==0
          y= x_GI;  % No channel
      else  % Multipath fading channel
        channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);
        h=zeros(1,Lch);
        h(Delay+1)=channel; % cir: channel impulse response
        y = conv(x_GI,h); 


%         [y, pathGains] = rayleighChan(x_GI');
%         chImp = randi([0 1],1000,1);
%         h = rayleighChan(chImp);
      end
      
      if i==0 % Only to measure the signal power for adding AWGN noise
        y1=y(1:Nframe*Nsym); 
        sigPow = sigPow + y1*y1';
        continue;
      end
      
      % Add AWGN noise________________________________________________
%       snr = EbN0(i)+10*log10(Nbps*(Nused/Nfft)); % SNR vs. Eb/N0
      snr = EbN0;
      noise_mag = sqrt((10.^(-snr/10))*sigPow/2);
      y_GI = y + noise_mag*(randn(size(y))+j*randn(size(y)));
      
      % Rx_____________________________________________________________
      kk1=(NgType==2)*Ng+[1:Nsym];
      kk2=1:Nfft;
      kk3=1:Nused;
      kk4=Nused/2+Nvc+1:Nfft;
      kk5=(Nvc~=0)+[1:Nused/2];
      
      if Ch==1
         H= fft([h zeros(1,Nfft-Lch)]); % Channel frequency response
         H_shift(kk3)= [H(kk4) H(kk5)]; 
      end
      
      for k=1:Nframe
         Y(kk2)= fft(remove_GI(Ng,Nsym,NgType,y_GI(kk1)));
         Y_shift=[Y(kk4) Y(kk5)];
         
         if Ch==0
             Xmod_r(kk3) = Y_shift;
         else
             Xmod_r(kk3)= Y_shift./H_shift;  % Equalizer - channel compensation
         end
         
         kk1=kk1+Nsym; 
         kk2=kk2+Nfft;
         kk3=kk3+Nused;
         kk4=kk4+Nfft; 
         kk5=kk5+Nfft;
      end
      
      X_r=qamdemod(Xmod_r*norms(Nbps),M,'gray');
      Neb=Neb+sum(sum(de2bi(X_r,Nbps)~=de2bi(X(:,1)',Nbps)));
      Ntb=Ntb+Nused*Nframe*Nbps;  %[Ber,Neb,Ntb]=ber(bit_Rx,bit,Nbps); 
      
      if Neb>Target_neb
          break;
      end
   end
   
   if i==0
     sigPow= sigPow/Nsym/Nframe/N_iter;
     fprintf('Signal power= %11.3e\n', sigPow);
%      fprintf(fid,'%%Signal power= %11.3e\n%%EbN0[dB]       BER\n', sigPow);
    else
     Ber = Neb/Ntb/10;     
     fprintf('EbN0=%3d[dB], BER=%4d/%8d =%11.3e , Throughput: %3d Mbps \n', EbN0(i), Neb,Ntb,Ber,throughput)
     
%      fprintf(fid, '%d\t%11.3e\n', EbN0(i), Ber);
     berVal(i) = Ber;
     if Ber<1e-6
         break;  
     end
   end
   
end
fprintf('SCS: %d, Cyclic Prefix: %d, Mod Order: %d, ', Nfft,Ng,Nbps);
   
% if (fid~=0)
%     fclose(fid);  
% end
% disp('Simulation is finished');
% plot_ber(file_name,Nbps);
end


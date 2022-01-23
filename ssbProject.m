clc;
clear all
close all

[audioTime,FS]=audioread('eric.wav');
%sound(audio,FS);
N = length(audioTime);       % no. of samples in audio signal audio

%% Specturm of the audio 

audioFreq = (fftshift(fft(audioTime))); 
F = linspace(-FS/2,FS/2,N);  %frequency domain of the spectrum 
plot(F,abs(audioFreq));
title('Spectrum of input audio')
xlim([-25000 25000]);

%%  ideal filter  

%first we should get the number of samples that lie between -4khz and 4khz
samples = 68542; % (N/48000) * 8000 (Total samples from -4khz to 4khz)
pass=ones(samples,1);
cutoff=zeros((N-samples)/2,1);
imp=[cutoff;pass;cutoff];
M = imp.*audioFreq; %filtered audio in freq. domain

figure
plot(F,abs(M));  %plot the spectrum of the filtered audio
xlim([-5000 5000]);
title('Spectrum of filtered audio')

m = real(ifft(ifftshift(M)));
t = linspace(0, length(m)/FS,length(m));
figure
plot(t,m);  %plot the filtered audio in time 
title('Filtered audio in time')
%sound(m,FS);

%% DSB-SC Modulation
FC = 100000;
FS_new = 5*FC;
m_resampled = resample(real(m),FS_new,FS);
t_resampled = linspace(0, length(m_resampled)/FS_new,length(m_resampled));
t_resampled = t_resampled' ; % transpose of time after resampling,
                             % to make the dimenstions the of resampled message
                             % and the carrier equal
carrier = cos(2*pi*FC*t_resampled);
u_sc = m_resampled.*carrier;
 
% frequency domain of modulated signal dsb-sc
U_sc = (fftshift(fft(u_sc))); 
fdomain_sc = linspace(-FS_new/2,FS_new/2,length(U_sc));


%SSB ideal filter get LSB (from -100k Hz to 100k Hz)
samples_LSB = 1713536 ;  %(N/48000)*200000
pass = ones(samples_LSB,1);  
cutoff = zeros((length(U_sc) - samples_LSB)/2,1);
imp = [cutoff;pass;cutoff];

ssbFilter = U_sc .* imp;
ssbFilterFreq = linspace(-FS_new/2,FS_new/2,length(ssbFilter));
figure
plot(ssbFilterFreq,abs(ssbFilter));
title('SSB-SC LSB  ideal low pass filter')
xlim([-150000 150000])

ssbFilterTime = real(ifft(ifftshift(ssbFilter)));
%% butter filter

Fn = FS_new/2;                             
wc1=(FC-4000)/Fn;
wc2=(FC)/Fn;
[b,a] = butter(4, [wc1 wc2],'bandpass');
filteredintime = filter(b,a,u_sc);
ssbButter=fftshift(fft(filteredintime));
ssbFilterFreq1 = linspace(-FS_new/2,FS_new/2,length(ssbButter));
plot(ssbFilterFreq1,abs(ssbButter));
title('SSB LSB butter-worth filter')

%% adding noise to SSB-SC signal
% snr  = 0 db
snr = 0;
ssbNoise0 = awgn(ssbFilterTime,snr,'measured');
% snr = 10 db
snr2 = 10;
ssbNoise10 = awgn(ssbFilterTime,snr2,'measured');
% snr = 30 db
snr3 = 30;
ssbNoise30 = awgn(ssbFilterTime,snr3,'measured');

CoherentDetector(FC,ssbNoise0,FS,FS_new,samples,U_sc,snr,0);
CoherentDetector(FC,ssbNoise10,FS,FS_new,samples,U_sc,snr2,0);
CoherentDetector(FC,ssbNoise30,FS,FS_new,samples,U_sc,snr3,0);

%% frequency error
CoherentDetector(100.1*1000,ssbNoise0,FS,FS_new,samples,U_sc,snr,0);
CoherentDetector(100.1*1000,ssbNoise10,FS,FS_new,samples,U_sc,snr2,0);
CoherentDetector(100.1*1000,ssbNoise30,FS,FS_new,samples,U_sc,snr3,0);

%% phase error = pi/9
 CoherentDetector(FC,ssbNoise0,FS,FS_new,samples,U_sc,snr,pi/9);
 CoherentDetector(FC,ssbNoise10,FS,FS_new,samples,U_sc,snr2,pi/9);
 CoherentDetector(FC,ssbNoise30,FS,FS_new,samples,U_sc,snr3,pi/9);


%% SSB-TC

% generating DSB-TC
Ac = 2*max(real(m));  %DC bias is twice the maximum of the message
u_tc = (Ac + m_resampled).*carrier;
U_tc = (fftshift(fft(u_tc)));
fdomain_tc = linspace(-FS_new/2,FS_new/2,length(U_tc));


pass = ones(samples_LSB  ,1);
cutoff = zeros((length(U_tc) - (samples_LSB))/2,1);
imp = [cutoff ;pass; cutoff];
ssbTCFilter = U_tc .* imp;
ssbTCFilterFreq = linspace(-FS_new/2,FS_new/2,length(ssbTCFilter));

figure
subplot(1,2,2);
plot(ssbTCFilterFreq,abs(ssbTCFilter));
xlim([0.95*10^5 1.05*10^5]);
ylim([0 4000]);
subplot(1,2,1)
plot(ssbTCFilterFreq,abs(ssbTCFilter));
xlim([-1.05*10^5 -0.95*10^5]);
ylim([0 4000]);

suptitle('SSB-TC LSB Ideal Low Pass Filter ')
ssbTCFilterTime=ifft(ifftshift(ssbTCFilter));

envelope_tc = abs(hilbert(ssbTCFilterTime)); 
envelope_tc = envelope_tc-abs(mean(envelope_tc)); % DC block
ssb_dmod_tc = resample(envelope_tc,FS,FS_new);    

t2 = linspace(0,length(ssb_dmod_tc)/FS,length(ssb_dmod_tc));

figure
plot(t2,abs(ssb_dmod_tc));
title('Recieved LSSB-TC in time');

%sound(ssb_dmod_tc,FS);



close all
clc
clear
[audio,FS]=audioread('eric.wav');

N = length(audio);       % no. of samples in audio signal audio

%% Specturm of the audio

audioFreq = (fftshift(fft(audio))); 
F = linspace(-FS/2,FS/2,N);  %frequency domain of the spectrum 
plot(F,abs(audioFreq));
title('Spectrum of input audio')
xlim([-25000 25000]);

%% ideal filter 

%first we should get the number of samples that lie between -4khz and 4khz
samples = 68542; % (N/48000) * 8000 (Total samples from -4khz to 4khz)
pass = ones(samples,1);
cutoff = zeros((N-samples)/2,1);
imp=[cutoff; pass; cutoff];

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
fdomain_sc = linspace(-FS_new/2,FS_new/2,length(u_sc));
figure
plot(fdomain_sc,abs(U_sc));
title('Spectrum of DSB-SC signal');


%% DSB-TC Modulation

Ac = 2*max(real(m));  %DC bias is twice the maximum of the message
u_tc = (Ac + m_resampled).*carrier;
U_tc = (fftshift(fft(u_tc)));
fdomain_tc = linspace(-FS_new/2,FS_new/2,length(u_tc));
figure
plot(fdomain_tc,abs(U_tc));
title(' Specturm of DSB-TC signal ')


%% Envelope detection w/o noise 
envelope_sc = abs(hilbert(u_sc));
envelope_tc1 = abs(hilbert(u_tc));

envelope_tc = envelope_tc1-abs(mean(envelope_tc1));

y1_dmod_sc = resample(envelope_sc,FS,FS_new);
y2_dmod_tc = resample(envelope_tc,FS,FS_new);

t1 = linspace(0,length(y1_dmod_sc)/FS,length(y1_dmod_sc));
t2 = linspace(0,length(y2_dmod_tc)/FS,length(y2_dmod_tc));

figure
plot(t1,abs(y1_dmod_sc));
title('Recieved DSB-SC in time');

figure
plot(t2,abs(y2_dmod_tc));
title('Recieved DSB-TC in time');

%sound(y1_dmod_sc,FS);
%sound(y2_dmod_tc,FS);

%observation: 
%the detection of the DSB-SC signal ouputs  a distorted audio
%the detection of the DSB-TC signal outputs normal audio (desired)
% the envelope detector can be used with DSB-TC



%% adding noise to DSB-SC signal

% snr  = 0 db
snr = 0;
dsbSCNoise0 = awgn(u_sc,snr,'measured');
% snr = 10 db
snr2 = 10;
dsbSCNoise10 = awgn(u_sc,snr2,'measured');
% snr = 30 db
snr3 = 30;
dsbSCNoise30 = awgn(u_sc,snr3,'measured');


%% DSB-SC with noise demodulation using coherent detection
% multiplying the modulated signal with same carrier, frequency error = 0
 CoherentDetector(FC,dsbSCNoise0,FS,FS_new,samples,u_sc,snr,0);
 CoherentDetector(FC,dsbSCNoise10,FS,FS_new,samples,u_sc,snr2,0);
 CoherentDetector(FC,dsbSCNoise30,FS,FS_new,samples,u_sc,snr3,0);

%% frequency error
 CoherentDetector(100.1*1000,dsbSCNoise0,FS,FS_new,samples,u_sc,snr,0);
 CoherentDetector(100.1*1000,dsbSCNoise10,FS,FS_new,samples,u_sc,snr2,0);
 CoherentDetector(100.1*1000,dsbSCNoise30,FS,FS_new,samples,u_sc,snr3,0); 
 
%% phase error = pi/9
 CoherentDetector(FC,dsbSCNoise0,FS,FS_new,samples,u_sc,snr,pi/9);
 CoherentDetector(FC,dsbSCNoise10,FS,FS_new,samples,u_sc,snr2,pi/9);
 CoherentDetector(FC,dsbSCNoise30,FS,FS_new,samples,u_sc,snr3,pi/9);





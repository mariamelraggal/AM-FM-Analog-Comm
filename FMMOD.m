clc
clear 
close all
%% read audio
[file,path] = uigetfile('*.wav');
[signal,fs]=audioread(file);
%sound(signal,fs);
n=length(signal);

%% calculate number of samples in LPF
samples=ceil((8000*n)/fs);
signalspectrum=fftshift(fft(signal));
signalfreqdomain=linspace(-fs/2,fs/2,n);
figure
plot(signalfreqdomain,signalspectrum);
title('Audio Signal spectrum');
x=ones(samples,1);
y=zeros(floor((n-samples)/2),1);
imp=[y;x;y];
filteredinfreq=imp.*signalspectrum;
figure
plot(signalfreqdomain,filteredinfreq);
xlim([-5000 5000]);
title('Filtered Signal spectrum');
filteredintime=real(ifft(ifftshift(filteredinfreq)));
timeFilterted = linspace(0,length(filteredintime)/fs,length(filteredintime));
figure
plot(filteredintime)
title('Filtered Signal in time');
mt=filteredintime;
%sound(filteredintime,fs);

%% modulation
fc=100*1000;
fs_new=5*fc;
mt=resample(mt,fs_new,fs);
Qt=max(cumsum(mt)*(1/(fs_new)));
mtt=cumsum(mt)*(1/(fs_new));
kf=floor(0.286/(Qt));
t=linspace(1,length(mt)/fs_new,length(mt));
t=t';
st=cos(2*pi*fc*t+2*pi*kf*mtt);
fmspectrum=fftshift(fft(st)/fs_new);
fmfreq=linspace(-fs_new/2,fs_new/2,length(fmspectrum));
figure
plot(fmfreq,fmspectrum);
title('Modulated FM signal spectrum');

%% demodulation
sigdiff=diff(st);
siged=0.000286*abs(hilbert(sigdiff*fs_new));
DCBlock=siged-mean(siged);
demsig=resample(DCBlock,fs,fs_new);
time = linspace(0,length(demsig)/fs,length(demsig));
figure
plot(time,demsig);
title('Demodulated FM signal in time');
ylim([-0.25 0.2])
%sound(demsig,fs);
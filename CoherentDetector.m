function [] = CoherentDetector(FC, message,FS,FS_new,samples,u_sc,snr,phase)
t = linspace(0,length(message)/FS_new,length(message));
carrier = 2*cos(2*pi*FC*t + phase);
carrier = carrier';
dsbSC_DEMOD = message.*carrier;
dsbSC_DEMOD_FREQ = fftshift(fft(dsbSC_DEMOD));

% applying LPF to extract the msg

pass=ones(samples,1);
cutoff =zeros((length(u_sc) - samples)/2,1);
imp=[cutoff;pass;cutoff];
M_DEMOD = imp.*dsbSC_DEMOD_FREQ;
f_demod_domain = linspace(-FS_new/2,FS_new/2,length(M_DEMOD));

figure
plot(f_demod_domain,M_DEMOD);
s = sprintf('Demoulated signal frequency domain snr = %d',snr);

title(s);
xlim([-5000 5000]);

dsbSC_demod_time = ifft(ifftshift(M_DEMOD));
msgDEMOD = resample(real(dsbSC_demod_time),FS,FS_new);
t3 = linspace(0,length(msgDEMOD)/FS, length(msgDEMOD));
figure
plot(t3,msgDEMOD);
s = sprintf('Demoulated signal time domain snr = %d',snr);
title(s);

%sound(msgDEMOD,FS);
end
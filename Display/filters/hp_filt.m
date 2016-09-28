function x = hp_filt(f_order, hf, sig, dt)

% this outputs the data after a high-pass filter has been used

sampling_rate = 1/dt;

Fnyquist = sampling_rate/2;

cutoff = hf/Fnyquist;

[b, a] = butter(f_order, cutoff); % butter filter

xfil = filtfilt(b, a, sig);

x = sig - xfil;











function x = lp_filt(f_order, lf, sig, dt)

% this outputs the data after a low-pass filter has been used

sampling_rate = 1/dt;

Fnyquist = sampling_rate/2;

cutoff = lf/Fnyquist;

[b, a] = butter(f_order, cutoff); % butter filter

x = filtfilt(b, a, sig);











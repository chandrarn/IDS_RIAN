function x = bp_filt(f_order, lf, hf, sig, dt)

% this outputs the data after a band-pass filter has been used

sampling_rate = 1/dt;

Fnyquist = sampling_rate/2;

lcutoff = lf/Fnyquist;
hcutoff = hf/Fnyquist;

[b1, a1] = butter(f_order, lcutoff); % butter filter
[b2, a2] = butter(f_order, hcutoff);

xfil1 = filtfilt(b1, a1, sig);
xfil2 = filtfilt(b2, a2, sig);

x = xfil2 - xfil1;











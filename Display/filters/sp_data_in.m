function [t, dt, x] = sp_data_in(sig)

% this for the raw data from the dtaq, not for the already processed data

x = mdsvalue(sig);
t_bad = mdsvalue(['dim_of(' sig ')']);
% the time base has an extra value so this needs to be done
t_length = length(x);

% the dt saved for the sp's is not correct ...
f_actual = 501253; %Hz
dt_actual = 1/f_actual;

[Y, I] = min(abs(t_bad));

t = [-1*((I - 1): -1: 1)*dt_actual + Y, Y, ...
    dt_actual*(1: (t_length - I)) + Y]';

dt = dt_actual;

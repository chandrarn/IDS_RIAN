function [t, dt, x] = gen_data_in(sig)

% this gets generic data from landau

x = mdsvalue(sig);
dt = mdsvalue(['samplinginterval(' sig ')']);
tmin = mdsvalue(['minval(dim_of(' sig '))']);
tlength = mdsvalue(['size(dim_of(' sig '))']);
t = tmin + dt*((1: tlength) - 1);


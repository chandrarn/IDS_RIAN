function [t, dt, x] = gen_data_in_smart(sig)

% this gets generic data from landau

x = mdsvalue(sig);
if isa(x, 'char')
    t = NaN;
    dt = NaN;
    x = NaN;
elseif isa(x, 'numeric')
    dt = mdsvalue(['samplinginterval(' sig ')']);
    tmin = mdsvalue(['minval(dim_of(' sig '))']);
    tlength = mdsvalue(['size(dim_of(' sig '))']);
    t = tmin + dt*((1: tlength) - 1);
end
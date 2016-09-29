function [t, dt, x] = gen_data_in(sig,HitTree)

% this gets generic data from landau

x = NATIVEvalue(HitTree.get(sig));
dt = NATIVEvalue(HitTree.get(['samplinginterval(' sig ')']));
tmin = NATIVEvalue(HitTree.get(['minval(dim_of(' sig '))']));
tlength = NATIVEvalue(HitTree.get(['size(dim_of(' sig '))']));

t = tmin + dt*(double(1: tlength) - 1);


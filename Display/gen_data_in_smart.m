function [t, dt, x] = gen_data_in_smart(sig,HitTree)

% this gets generic data from landau

x = mdsvalue(sig);
if isa(x, 'char')
    t = NaN;
    dt = NaN;
    x = NaN;
elseif isa(x, 'numeric')
    dt = NATIVEvalue(HitTree.get(['samplinginterval(' sig ')']));
    tmin = NATIVEvalue(HitTree.get(['minval(dim_of(' sig '))']));
    tlength = NATIVEvalue(HitTree.get(['size(dim_of(' sig '))']));
    t = tmin + dt*((1: tlength) - 1);
end
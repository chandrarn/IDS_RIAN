function [cb_field, ins_depth] = imp_data_in_raw(array, ...
    tin)

% this gets raw IMP data from landau

probe = ['01'; '02'; '03'; '04'; '05'; '06'; ...
    '07'; '08'; '09'; '10'; '11'; '12'; ...
    '13'; '14'; '15'; '16'; '17'];

N = length(probe(:,1));

dir = ['R'; 'P'; 'T'];
% r is in the i direction
% p is in the j direction
% t is in the k direction

cb_field = zeros(length(dir), N, length(tin));
ins_depth = zeros(1, N);

for j = 1: N
    for jj = 1: length(dir)
        if jj == 1 % the radial probes don't have a b_winding
            node = ['\B_IMP_' array '_' dir(jj, :) probe(j, :)];
        elseif jj >= 2
            node = ['\B_IMP_' array '_' dir(jj, :) probe(j, :) ':b_winding'];
        end
        dt = mdsvalue(['samplinginterval(' node ')']);
        tmin = mdsvalue(['minval(dim_of(' node '))']);
        tlength = mdsvalue(['size(dim_of(' node '))']);
        if isa(tmin, 'char')
            cb_field(jj,j,:) = zeros(1, length(tin))*NaN;
        elseif isa(tmin, 'numeric')
            b_field = mdsvalue(node);
            t = tmin + dt*((1: tlength) - 1);
            cb_field(jj,j,:) = interp1(t, b_field, tin);
        end
    end
    ins_depth(j) = mdsvalue(['\B_IMP_M_R' probe(j, :) 'R']);
end

function [cb_field, ins_depth] = imp_data_in_construct(array, ...
    tin, filterdata, f_order, low_filt, high_filt)

% using a slanted baseline and matrix

% this gets the IMP data from landau

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
    pnode = ['\B_IMP_' array '_P' probe(j, :) ':b_winding'];
    tnode = ['\B_IMP_' array '_T' probe(j, :) ':b_winding'];
    rnode = ['\B_IMP_' array '_R' probe(j, :)];
    dtp = mdsvalue(['samplinginterval(' pnode ')']);
    dtt = mdsvalue(['samplinginterval(' tnode ')']);
    dtr = mdsvalue(['samplinginterval(' rnode ')']);
    tminp = mdsvalue(['minval(dim_of(' pnode '))']);
    tmint = mdsvalue(['minval(dim_of(' tnode '))']);
    tminr = mdsvalue(['minval(dim_of(' rnode '))']);
    tlengthp = mdsvalue(['size(dim_of(' pnode '))']);
    tlengtht = mdsvalue(['size(dim_of(' tnode '))']);
    tlengthr = mdsvalue(['size(dim_of(' rnode '))']);
    if isa(tminr, 'char')
        cbw_rad = zeros(1, length(tin))*NaN;
    elseif isa(tminr, 'numeric')
        trad = tminr + dtr*((1: tlengthr) - 1);
        bw_rad = mdsvalue(['slanted_baseline2(sub_baseline_string("\' rnode '"))']);
        if filterdata == 2;
            % bandpass filter
            bw_fil_r = bp_filt(f_order, low_filt, ...
                high_filt, bw_rad, dtr);
        elseif filterdata == 1;
            % lowpass filter
            bw_fil_r = lp_filt(f_order, low_filt, ...
                bw_rad, dtr);
        elseif filterdata == 0;
            bw_fil_r = bw_rad;
        elseif filterdata == 3;
            % filter at half Nyquist
            bw_fil_r = lp_filt(f_order, 1/dtr/4, ...
                bw_rad, dtr);
        elseif filterdata == 4;
            % inverse bandpass filter
            bw_fil_r = ibp_filt(f_order, low_filt, ...
                high_filt, bw_rad, dtr);
        end
        cbw_rad = interp1(trad, bw_fil_r, tin);
    end
    if isa(tminp, 'char')
        cbw_pol = zeros(1, length(tin))*NaN;
    elseif isa(tminp, 'numeric')
        tpol = tminp + dtp*((1: tlengthp) - 1);
        bw_pol = mdsvalue(['slanted_baseline2(sub_baseline_string("\' pnode '"))']);
        if filterdata == 2;
            % bandpass filter
            bw_fil_p = bp_filt(f_order, low_filt, ...
                high_filt, bw_pol, dtp);
        elseif filterdata == 1;
            % lowpass filter
            bw_fil_p = lp_filt(f_order, low_filt, ...
                bw_pol, dtp);
        elseif filterdata == 0;
            bw_fil_p = bw_pol;
        elseif filterdata == 3;
            % filter at half Nyquist
            bw_fil_p = lp_filt(f_order, 1/dtp/4, ...
                bw_pol, dtp);
        elseif filterdata == 4;
            % inverse bandpass filter
            bw_fil_p = ibp_filt(f_order, low_filt, ...
                high_filt, bw_pol, dtp);
        end
        cbw_pol = interp1(tpol, bw_fil_p, tin);
    end
    if isa(tmint, 'char')
        cbw_tor = zeros(1, length(tin))*NaN;
    elseif isa(tmint, 'numeric')
        ttor = tmint + dtt*((1: tlengtht) - 1);
        bw_tor = mdsvalue(['slanted_baseline2(sub_baseline_string("\' tnode '"))']);
        if filterdata == 2;
            % bandpass filter
            bw_fil_t = bp_filt(f_order, low_filt, ...
                high_filt, bw_tor, dtt);
        elseif filterdata == 1;
            % lowpass filter
            bw_fil_t = lp_filt(f_order, low_filt, ...
                bw_tor, dtt);
        elseif filterdata == 0;
            bw_fil_t = bw_tor;
        elseif filterdata == 3;
            % filter at half Nyquist
            bw_fil_t = lp_filt(f_order, 1/dtt/4, ...
                bw_tor, dtt);
        elseif filterdata == 4;
            % inverse bandpass filter
            bw_fil_t = ibp_filt(f_order, low_filt, ...
                high_filt, bw_tor, dtt);
        end
        cbw_tor = interp1(ttor, bw_fil_t, tin);
    end
    cb_field(1,j,:) = cbw_rad;
    if j == 1 % there is no rot ang for the 1st probe
        % b/c there is no toroidal probe
        cb_field(2,j,:) = cbw_pol;
        cb_field(3,j,:) = cbw_tor;
    elseif j > 1
        rotnode = ['\B_IMP_' array '_T' probe(j, :) ':ROT_ANG'];
        rot_ang = mdsvalue(rotnode);
        cb_field(2,j,:) = cbw_pol*cos(rot_ang) - cbw_tor*sin(rot_ang);
        cb_field(3,j,:) = cbw_tor*cos(rot_ang) + cbw_pol*sin(rot_ang);
    end
    ins_depth(j) = mdsvalue(['\B_IMP_M_R' probe(j, :) 'R']);
end

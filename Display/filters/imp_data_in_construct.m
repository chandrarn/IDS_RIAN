function [cb_field, ins_depth] = imp_data_in_construct(array, ...
    tin, filterdata, f_order, low_filt, high_filt, shot, shift, dafi)

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

% calibration factors based on dafi impedence.
% these values good from shot 121973
pol_scale = [1.02;
    1.016;
    1.014;
    1.018;
    1.013;
    1.022;
    1.027;
    1.024;
    1.025;
    1.022;
    1.02;
    1.028;
    1.029;
    1.022;
    1.022;
    1.024;
    1.021];

% there is no first toroidal probe
tor_scale = [0;
    1.025;
    1.012;
    1.014;
    1.008;
    1.022;
    1.024;
    1.024;
    1.023;
    1.012;
    1.019;
    1.02;
    1.022;
    1.019;
    1.021;
    1.026;
    1.025];


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
        % shifting time base here for digi differences
        if strcmp(shift, 'shift_time') == 1
            trad = imp_time_shift(trad, shot, array, probe(j,:), 'R');
        end
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
        % shifting time base here for digi differences
        if strcmp(shift, 'shift_time') == 1
            tpol = imp_time_shift(tpol, shot, array, probe(j,:), 'P');
        end
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
        if strcmp(dafi, 'dafi_cf') == 1
            if shot >= 121973
                cbw_pol = cbw_pol*pol_scale(j);
            end
        end
    end
    if isa(tmint, 'char')
        cbw_tor = zeros(1, length(tin))*NaN;
    elseif isa(tmint, 'numeric')
        ttor = tmint + dtt*((1: tlengtht) - 1);
        % shifting time base here for digi differences
        if strcmp(shift, 'shift_time') == 1
            ttor = imp_time_shift(ttor, shot, array, probe(j,:), 'T');
        end
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
        if strcmp(dafi, 'dafi_cf') == 1
            if shot >= 121973
                cbw_tor = cbw_tor*tor_scale(j);
            end
        end
    end
    cb_field(1,j,:) = cbw_rad;
    if strcmp(array, 'M') == 1
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
    else
        rotnode = ['\B_IMP_' array '_T' probe(j, :) ':ROT_ANG'];
        rot_ang = mdsvalue(rotnode);
        if isa(rot_ang, 'char')
            cb_field(2,j,:) = cbw_pol;
            cb_field(3,j,:) = cbw_tor;
        elseif isa(tmint, 'numeric')
            cb_field(2,j,:) = cbw_pol*cos(rot_ang) - cbw_tor*sin(rot_ang);
            cb_field(3,j,:) = cbw_tor*cos(rot_ang) + cbw_pol*sin(rot_ang);
        end
    end
    ins_depth(j) = mdsvalue(['\B_IMP_M_R' probe(j, :) 'R']);
end

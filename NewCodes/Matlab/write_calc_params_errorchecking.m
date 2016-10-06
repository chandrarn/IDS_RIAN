% calculate parameters and write them to the analysis tree for HIT-SI3
%5/10/2016 - Added error handling in case a signal fails to store - KDM
function write_calc_params_errorchecking(shot)

[a, status] = mdsopen('landau.hit::hitsi3', shot);

do_Iq = 1;
do_Vq = 1;
do_Pq = 1;

try
[t_iquad, dt_iquad, sig_iquad] = ...
    gen_data_in(['sqrt( sigadd( sigadd( sigmul(\i_inj_a, \i_inj_a), ' ...
    'sigmul(\i_inj_b, \i_inj_b) ), sigmul(\i_inj_c, \i_inj_c) ) )']);
catch
do_Iq=0;
fprintf('Unable to get current values\n');
end

try
[t_vquad, dt_vquad, sig_vquad] = ...
    gen_data_in(['sqrt( sigadd( sigadd( sigmul(\v_inj_a, \v_inj_a), ' ...
    'sigmul(\v_inj_b, \v_inj_b) ), sigmul(\v_inj_c, \v_inj_c) ) )']);
catch
do_Vq = 0;
fprintf('Unable to get voltage values\n');
end

try
[t_psiquad, dt_psiquad, sig_psiquad] = ...
    gen_data_in(['sqrt( sigadd( sigadd( sigmul(\psi_inj_a, \psi_inj_a), ' ...
    'sigmul(\psi_inj_b, \psi_inj_b) ), sigmul(\psi_inj_c, \psi_inj_c) ) )']);
catch
do_Pq=0;
fprintf('Unable to get flux values\n');
end
inj_names = ['a', 'b', 'c'];

for j = 1: length(inj_names)
    
    if do_Vq && do_Pq
    [t_kdot.(inj_names(j)), dt_kdot.(inj_names(j)), sig_kdot.(inj_names(j))] = ...
        gen_data_in(['2 * sigmul(\psi_inj_' inj_names(j) ',\v_inj_' inj_names(j) ')']);
    end
    
    if do_Vq && do_Iq
    [t_p.(inj_names(j)), dt_p.(inj_names(j)), sig_p.(inj_names(j))] = ...
        gen_data_in(['sigmul(\i_inj_' inj_names(j) ',\v_inj_' inj_names(j) ')']);
    end
    
    if do_Iq && do_Pq
    [t_lambda.(inj_names(j)), dt_lambda.(inj_names(j)), sig_lambda.(inj_names(j))] = ...
        gen_data_in(['4 * $pi * 1e-7 * sqrt(sigdiv(sihi_smooth(sigmul(\i_inj_' inj_names(j) ...
        ',\i_inj_' inj_names(j) ')), sihi_smooth(sigmul(\psi_inj_' inj_names(j) ',\psi_inj_' ...
        inj_names(j) '))))']);
    end
    
end

%------------------------------------------------------------------------%
%------------------------------------------------------------------------%

% Write data to tree

if do_Iq
mdsput('\I_INJ_QUAD', 'build_signal(build_with_units($1,"A"), *, build_with_units($2,"s"))', ...
    sig_iquad, t_iquad);
end
if do_Vq
mdsput('\V_INJ_QUAD', 'build_signal(build_with_units($1,"V"), *, build_with_units($2,"s"))', ...
    sig_vquad, t_vquad);
end
if do_Pq
mdsput('\PSI_INJ_QUAD', 'build_signal(build_with_units($1,"Wb"), *, build_with_units($2,"s"))', ...
    sig_psiquad, t_psiquad);
end

for j = 1: length(inj_names)
    
    if do_Vq && do_Pq
    mdsput(['\KDOT_INJ_' inj_names(j)], ...
        'build_signal(build_with_units($1,"Wb^2"), *, build_with_units($2,"s"))', ...
        sig_kdot.(inj_names(j)), t_kdot.(inj_names(j)));
    end
    
    if do_Vq && do_Iq
    mdsput(['\P_INJ_' inj_names(j)], ...
        'build_signal(build_with_units($1,"W"), *, build_with_units($2,"s"))', ...
        sig_p.(inj_names(j)), t_p.(inj_names(j)));
    end
    
    if do_Iq && do_Pq
    mdsput(['\LAMBDA_INJ_' inj_names(j)], ...
        'build_signal(build_with_units($1,"m^-1"), *, build_with_units($2,"s"))', ...
        sig_lambda.(inj_names(j)), t_lambda.(inj_names(j)));
    end
end

% signals made from signals made earlier in this code
if do_Vq && do_Pq
[t_kdot_tot, dt_kdot_tot, sig_kdot_tot] = ...
    gen_data_in('sigadd(sigadd(\KDOT_INJ_A, \KDOT_INJ_B), \KDOT_INJ_C)');
end

if do_Vq && do_Iq
[t_p_tot, dt_p_tot, sig_p_tot] = ...
    gen_data_in('sigadd(sigadd(\P_INJ_A, \P_INJ_B), \P_INJ_C)');
end

if do_Iq && do_Pq
[t_lambda_quad, dt_lambda_quad, sig_lambda_quad] = ...
    gen_data_in('4 * $pi * 1e-7 * sigdiv(\i_inj_quad,\psi_inj_quad)');
end
%------------------------------------------------------------------------%
%------------------------------------------------------------------------%

% Write data to tree
if do_Vq && do_Pq
mdsput('\KDOT_INJ', 'build_signal(build_with_units($1,"Wb^2"), *, build_with_units($2,"s"))', ...
    sig_kdot_tot, t_kdot_tot);
end
if do_Vq && do_Iq
mdsput('\P_INJ', 'build_signal(build_with_units($1,"W"), *, build_with_units($2,"s"))', ...
    sig_p_tot, t_p_tot);
end
if do_Iq && do_Pq
mdsput('\LAMBDA_QUAD', 'build_signal(build_with_units($1,"m^-1"), *, build_with_units($2,"s"))', ...
    sig_lambda_quad, t_lambda_quad);
end

mdsclose;

%------------------------------------------------------------------------%
%------------------------------------------------------------------------%

end
%% bring the data in with a correct time base
function [t, dt, x] = gen_data_in(sig)

% this gets generic data from landau
x = mdsvalue(sig);
dt = mdsvalue(['samplinginterval(' sig ')']);
tmin = mdsvalue(['minval(dim_of(' sig '))']);
tlength = mdsvalue(['size(dim_of(' sig '))']);
tlength = cast(tlength,'like',dt);
t = tmin + dt*((1: tlength) - 1);

end


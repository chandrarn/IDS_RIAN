% Plot Data, mostly used now to process IDS data for later viewing by
% NIMcompare_1 or similar program.  Modified to incorporate frame adding in
% a spatial or temporal pattern.


% NB: Important Input settings Below %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

shots = [12981010, 12981710, 12981910, 12982010, 1294992, 1294993, 1294994]; % shot numbers
shotRef = [129810, 129817, 129819, 129820, 129499, 129499, 129499]; % if 'shots' is a false shot number, use this to pull cal data
                  % Set to 'NaN' if not using false shot number
line = [3]; % spectral line number, from longest to shortest wavelength

label1 = 'Impact Parameter [cm]'; % Toroidal Section
% label1 = 'Major Radius [cm]'; % Poloidal Section, fibers facing each other
label2 = 'Major Radius [cm]'; % Poloidal Section
% label1 = 'Displacement [cm]'; % Orthogonal Section

frame_sum = 1; % number of frames to sum (typically 1)
t_override = 0; % If 't_override' is set to 0, the time
% to analyze will be determined by 'shotDuration.m', or 't_override' can be a
% two value vector with start and end times in ms.
plots = 0; % makes plots for debugging

s.spSum = 0; % add chords together for better statistics
%     s.chSum = [10:14; 23:27]; % channels to add together.  The columns of each row are added together.
s.chSum = [9:11; 12:14; 16:18; 19:21; 23:25; 26:28];

s.tSum = 0; % add frames together in a temporal pattern
s.sumPatt = [50, 12, 200]; % [time point to start from, length of pattern, end]

addpath('S:\MC_IDS\Matlab Code\Core Codes');
addpath('S:\General Matlab');

for n = 1:length(shots)
    tic
    [out, param] = load_shot(shots(n), shotRef(n), line, frame_sum, t_override, plots, s);
    t = toc;
    disp(['Elapsed time is ' num2str(floor(t/60)) ' minutes, ' num2str(rem(t,60),2) ' seconds']);
    disp(['saving shot ' num2str(shots(n))]);
    
    mdsopen('landau.hit::hitsi', param.shotRef);
    [iinjxTime, dummy, iinjx] = gen_data_in('\I_INJ_X');
    dat.iinjxTime = 1e3 * iinjxTime;
    dat.iinjx = 1e-3 * iinjx;
    [iinjyTime, dummy, iinjy] = gen_data_in('\I_INJ_Y');
    dat.iinjyTime = 1e3 * iinjyTime;
    dat.iinjy = 1e-3 * iinjy;
    [ItorTime, dummy, Itor] = gen_data_in('\I_TOR_SPAAVG');
    dat.ItorTime = 1e3 * ItorTime;
    dat.Itor = 1e-3 * Itor;
    dat.impacts = mdsvalue('\IDS_IMPACTS');
    mdsclose();
    
    dat.time = out.time;
    dat.temp = out.temp;
    dat.vel = 1e-3*out.vel; % save in km/s
    dat.int = out.int;
    dat.peaks = param.peaks(:, 1);
    dat.raw = out.data; % raw data points
    dat.fit_par = out.fit_par; % Curve fit parameters [time, channel, 1:6]
    dat.param = param; % param structure
    dat.bounds = out.bounds; % Bounds for curve fitting grid, (time) x (channel) x (x1 x2 y1 y2)
    dat.guesses = out.guesses; % initial guess parameters, same size as fit_par
    dat.shotRef = param.shotRef; % real shot number
    dat.title = ['Shot ' num2str(dat.shotRef)];
    dat.label1 = label1;
    dat.label2 = label2;
    if param.calcError
        dat.dparam = out.dparams; % parameter uncertainty structure
        dat.stddev = out.stddev; % standard deviation of every data point
        dat.tempU = out.tempU; % upper error bar for temperature
        dat.tempL = out.tempL;
        dat.velU = 1e-3*out.velU; % upper error bar for velocity
        dat.velL = 1e-3*out.velL;
    end
    if exist('out.resnorm', 'var') % DNE if using LM curve fitting method
        dat.residual = sqrt(out.resnorm); % residual
    end
    if exist('out.exp', 'var')
        dat.exp = out.exp;
        dat.corrected = out.corrected;
    end
    cd('S:\MC_IDS\Matlab Code\Data Repository');
    if and(s.spSum, ~s.tSum)
        dat.chSum = s.chSum;
        save(['dat' num2str(shots(n)) 's'], 'dat');
    elseif and(~s.spSum, s.tSum)
        save(['dat' num2str(shots(n)) 't'], 'dat');
    elseif and(s.spSum, s.tSum)
        dat.chSum = s.chSum;
        save(['dat' num2str(shots(n)) 'st'], 'dat');
    else
        save(['dat' num2str(shots(n))], 'dat') % save 'dat' structure as 'dat<shot>.mat'
    end
end
    
    
    
    
    
    

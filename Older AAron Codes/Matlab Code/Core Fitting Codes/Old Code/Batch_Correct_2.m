% Plot Data, mostly used now to process IDS data for later viewing by
% NIMcompare_1 or similar program.  Modified to incorporate frame adding in
% a spatial or temporal pattern.


% NB: Important Input settings Below %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

import MDSplus.*
mdsconnect('landau.hit')


addpath('~/IDS/Matlab/');
try
    addAllThePaths;
end

shots = [129810]; % shot numbers
nBDmodes = [10]; % number of modes to save after BD filtering. Leave blank for no BD filtering
timeBound = [70:300]; % time point bounds for BD filtering.  Leave blank to use whole movie
line = [3]; % spectral line number, from longest to shortest wavelength

%Is the file/tree stucture from HIT-SI, or HIT-SI3?
hitsi3=0;
if (length(shots(1))>=6) 
    hitsi3=1;
end

    
label1 = 'Impact Parameter [cm]'; % Toroidal Section
% label1 = 'Major Radius [cm]'; % Poloidal Section, fibers facing each other
label2 = 'Major Radius [cm]'; % Poloidal Section
% label1 = 'Displacement [cm]'; % Orthogonal Section

frame_sum = 1; % number of frames to sum (typically 1)
plots = 0; % makes plots for debugging

s.spSum = 0; % add chords together for better statistics
%     s.chSum = [10:14; 23:27]; % channels to add together.  The columns of each row are added together.
s.chSum = [9:11; 12:14; 16:18; 19:21; 23:25; 26:28];

s.tSum = 0; % add frames together in a temporal pattern
s.sumPatt = [50, 12, 200]; % [time point to start from, length of pattern, end]

addpath('T:\RChandra\A-A-Ron Code\Matlab Code\Core Fitting Codes');
addpath('T:\RChandra\A-A-Ron Code\General Matlab');

for n = 1:length(shots)
    tic
    [out, param] = load_shot(shots(n), nBDmodes, timeBound, line, frame_sum, plots, s,hitsi3);
    t = toc;
    disp(['Elapsed time is ' num2str(floor(t/60)) ' minutes, ' num2str(rem(t,60),2) ' seconds']);
    disp(['saving shot ' num2str(shots(n))]);
    
    
    if(~hitsi3)
        %I ASSUME THAT THIS IS HOW TAGS WORK
        HitTree = Tree('hitsi',shots(n));
        Data = HitTree.getNode('\\I_INJ_X');
        [iinjxTime, dummy, iinjx] = Data.getData();
        dat.iinjxTime = 1e3 * iinjxTime;
        dat.iinjx = 1e-3 * iinjx;
        Data = HitTree.getNode('\\I_INJ_Y');
        [iinjyTime, dummy, iinjy] = Data.getData();
        dat.iinjyTime = 1e3 * iinjyTime;
        dat.iinjy = 1e-3 * iinjy;
        Data = HitTree.getNode('\\I_TOR_SPAAVG');
        [ItorTime, dummy, Itor]=Data.getData();
        dat.ItorTime = 1e3 * ItorTime;
        dat.Itor = 1e-3 * Itor;
        Data = HitTree.getNode('\\IDS_IMPACTS');
        dat.impacts=Data.getData();
        mdsclose();
%         % Also try getting data directly, with tags?
%         %mdsopen('landau.hit::hitsi', param.shotRef);
%         [iinjxTime, dummy, iinjx] = gen_data_in('\I_INJ_X');
%         dat.iinjxTime = 1e3 * iinjxTime;
%         dat.iinjx = 1e-3 * iinjx;
%         [iinjyTime, dummy, iinjy] = gen_data_in('\I_INJ_Y');
%         dat.iinjyTime = 1e3 * iinjyTime;
%         dat.iinjy = 1e-3 * iinjy;
%         [ItorTime, dummy, Itor] = gen_data_in('\I_TOR_SPAAVG');
%         dat.ItorTime = 1e3 * ItorTime;
%         dat.Itor = 1e-3 * Itor;
%         dat.impacts = mdsvalue('\IDS_IMPACTS');
%         mdsclose();
    else 
        HitTree = Tree('hitsi3',shot);
        Data = HitTree.getNode('\\I_INJ_A');
        [iinjaTime, dummy, iinja] = Data.getData();
        dat.iinjaTime = 1e3 * iinjaTime;
        dat.iinja = 1e-3 * iinja;
        Data = HitTree.getNode('\\I_INJ_B');
        [iinjbTime, dummy, iinjb] = Data.getData();
        dat.iinjbTime = 1e3 * iinjbTime;
        dat.iinjb = 1e-3 * iinjb;
         Data = HitTree.getNode('\\I_INJ_C');
        [iinjcTime, dummy, iinjc] = Data.getData();
        dat.iinjcTime = 1e3 * iinjcTime;
        dat.iinjc = 1e-3 * iinjc;
        Data = HitTree.getNode('\\I_TOR_SPAAVG');
        dat.ItorTime = 1e3 * ItorTime;
        dat.Itor = 1e-3 * Itor;
        Data = HitTree.getNode('\\IDS_IMPACTS');
        dat.impacts=Data.getData();
        mdsclose();
    end
    
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
    
    %Average chord velocity by injector phase. WHERE TO SAVE?
    [Phase,Velocity]=VelocityPhase(dat,[1,1.8],6)
    cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data');
    save(['PhaseTemp' num2str(shots(n)) ], '[Phase,Velocity]');
    
    cd('T:\IDS\Data Repository'); % IS THIS RIGHT?
    if and(s.spSum, ~s.tSum)
        dat.chSum = s.chSum;
        save(['dat' num2str(shots(n)) 's'], 'dat');
    elseif and(~s.spSum, s.tSum)
        save(['dat' num2str(shots(n)) 't'], 'dat');
    elseif and(s.spSum, s.tSum)
        dat.chSum = s.chSum;
        save(['dat' num2str(shots(n)) 'st'], 'dat');
    else% changed to temp: want to check difference between NSTX cine and mine
        save(['TEMP' num2str(shots(n)) num2str(nBDmodes)], 'dat') % save 'dat' structure as 'dat<shot><nBDmodes>.mat'
    end
end
    
    
    
    
    
    

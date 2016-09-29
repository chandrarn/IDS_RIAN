% Plot Data, mostly used now to process IDS data for later viewing by
% NIMcompare_1 or similar program.  Modified to incorporate frame adding in
% a spatial or temporal pattern. Aaron Hossack
% UPDATED TO OBJECT ORIENTED MDSplus BY RIAN CHANDRA, WNTR 2014
% UPDATED TO NEW CINE TYPE BY RIAN CHANDRA, SPR 2014

% NOTE: When using 129499 poloidal, need to manually put an extra zero at
% the end durring load shot, because 1294990 doesnt exist in the tree
% NB: Important Input settings Below %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
profile on;
import MDSplus.*
%mdsconnect('landau.hit')

try addpath('~/IDS/Matlab/');end
try
    addAllThePaths;
end

%129810,129819,129817,129820
%140708018,
%129810: 14.5
%128585: 53.5
%129213: 68.5
%shots = [150310021:150310023];%,
%shots = [150310023]; %,129530,129810 shot numbers
shots = [150601012:150601019];
nBDmodes = [10]; % number of modes to save after BD filtering. Leave blank for no BD filtering
timeBound = []; % time point bounds for BD filtering.  Leave blank to use whole movie
line = [3]; % spectral line number, from longest to shortest wavelength
useTree = 1; % Don't use ONLY in cases where caliibration was unable to save to the tree
%Is the file/tree stucture from HIT-SI, or HIT-SI3?
% (length(shots(1))>6)
hitsi3 = 0;
if (length(num2str(shots(1)))>8) 
    hitsi3 = 1;
else
    hitsi3 = 0;
end
%automatically trim the data file: ( leave blank if no
auto=[]; % time in mS for filtering. Leave blank for whole thing
%[how many points have to be above threshold, what the threshold is,
%[5,90,20]             % how many frames they have to be above it for

simulation = []; % leave emptpty unless using the impacts from an existing shot for simulation.
label1 = 'Impact Parameter [cm]'; % Toroidal Section
% label1 = 'Major Radius [cm]'; % Poloidal Section, fibers facing each other
label2 = 'Major Radius [cm]'; % Poloidal Section
% label1 = 'Displacement [cm]'; % Orthogonal Section

frame_sum = 1; % number of frames to sum (typically 1)
plots = 0; % makes plots for debugging

s.sim = []; % tag indicating type of simulation
            % 5 = PSI-TET, n^2 weighting and temperature
            % 7 = PSI-TET, n^2 weighting and temperature
            % 8 = NIMROD, n^2 weighting and temperature
           % LEAVE EMPTY [] IF NOT DOING SIMULATION
simDat = 'G:\Git\Visit\post_process_output.mat';
%simDat = 'C:\Documents and Settings\chandrarn\Desktop\Git\Visit\24x24_14_injRamp.mat';
%simDat = 'C:\Documents and Settings\chandrarn\Desktop\Git\Visit\PSITET_14kHz_xMHD.mat';
%simDat = 'G:\PSITET_14kHz_xMHD.mat';
% data file holding the currents/etc for the simulation

s.spSum = 0; % add chords together for better statistics
%     s.chSum = [10:14; 23:27]; % channels to add together.  The columns of each row are added together.
s.chSum = [9:11; 12:14; 16:18; 19:21; 23:25; 26:28];

s.tSum = 0; % add frames together in a temporal pattern
s.sumPatt = [50, 12, 200]; % [time point to start from, length of pattern, end]

addpath('T:\RChandra\A-A-Ron Code\Matlab Code\Core Fitting Codes');
addpath('T:\RChandra\A-A-Ron Code\General Matlab');
addpath('T:\IDS\Data Analysis');

for n = 1:length(shots)
    %try % DO ALL THE SHOT!!!!!!!!
    tic
%     if hitsi3
%         HitTree = Tree('analysis3',shots(n));
%     else
%         HitTree = Tree('analysis',shots(n));
%     end
%     Data = HitTree.getNode('\IDS_IMPACTS');
%     NATIVEvalue(Data.getData().data);
    cd('T:\IDS\Data Analysis');
    [out, param] = load_shot(shots(n), nBDmodes, timeBound, line, frame_sum, plots, s,useTree,auto);
    t = toc;
    disp(['Elapsed time is ' num2str(floor(t/60)) ' minutes, ' num2str(rem(t,60),2) ' seconds']);
    if ~isempty(s.sim)        
        disp(['saving shot ' num2str(s.sim(n)) num2str(shots(n))]);
    else
        disp(['saving shot ' num2str(shots(n))]);
    end
    
    if(~hitsi3)
        if isempty(s.sim)
            HitTree = Tree('hitsi',shots(n));
            Data = HitTree.getNode('\I_INJ_X');
            iinjx=Data.getData().data.getDoubleArray;
            %[iinjxTime, dummy, iinjx] = Data.getData();
            iinjxTime =  NATIVEvalue(Data.getDimensionAt(0).data());
            dat.iinjxTime = 1e3 * iinjxTime;
            dat.iinjx = 1e-3 * iinjx;
            Data = HitTree.getNode('\I_INJ_Y');
            iinjy=Data.getData().data.getDoubleArray;
            iinjyTime =  NATIVEvalue(Data.getDimensionAt(0).data());
            %[iinjyTime, dummy, iinjy] = Data.getData();
            dat.iinjyTime = 1e3 * iinjyTime;
            dat.iinjy = 1e-3 * iinjy;
            Data = HitTree.getNode('\I_TOR_SPAAVG');
            %[ItorTime, dummy, Itor]=Data.getData();
            Itor=Data.getData().data.getDoubleArray;
            ItorTime =  NATIVEvalue(Data.getDimensionAt(0).data());
            dat.ItorTime = 1e3 * ItorTime;
            dat.Itor = 1e-3 * Itor;

            mdsclose();
        else
            load(simDat);
            try % nimrod
                dat.iinjxTime = nimsave.time.*1e3;
                dat.iinjyTime = nimsave.time.*1e3;
                dat.ItorTime = nimsave.time.*1e3;
                dat.iinjx = nimsave.inj.I1.*1e-3;
                dat.iinjy = nimsave.inj.I2.*1e-3;
                dat.Itor = nimsave.Bprob.surf.Itor_spa_avg.*1e-3;
            catch % psitet
                dat.Itor = psisave.itor;
                dat.ItorTime = psisave.time;
            end
        end
            
        
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
        if isempty(s.sim)
            HitTree = Tree('hitsi3',shots(n));
            Data = HitTree.getNode('\I_INJ_A');
            dat.iinja = 1e-3 * NATIVEvalue(Data.getData.data());
            dat.iinjaTime = 1e3 *  NATIVEvalue(Data.getDimensionAt(0).data());
            Data = HitTree.getNode('\I_INJ_B');
            dat.iinjb = 1e-3 *NATIVEvalue(Data.getData.data());
            dat.iinjbTime = 1e3 *  NATIVEvalue(Data.getDimensionAt(0).data());
            Data = HitTree.getNode('\I_INJ_C');
            dat.iinjc = 1e-3 *NATIVEvalue(Data.getData.data());
            dat.iinjcTime = 1e3 *  NATIVEvalue(Data.getDimensionAt(0).data());
            Data = HitTree.getNode('\I_TOR_SPAAVG');
            Itor=Data.getData().data.getDoubleArray;
            ItorTime =  NATIVEvalue(Data.getDimensionAt(0).data());
            dat.ItorTime = 1e3 * ItorTime;
            dat.Itor = 1e-3 * Itor;

            mdsclose();
        else
            load(simDat);
            try % nimrod
                dat.iinjaTime = nimsave.time.*1e3;
                dat.iinjbTime = nimsave.time.*1e3;
                dat.iinjcTime = nimsave.time.*1e3;
                dat.ItorTime = nimsave.time.*1e3;
                dat.iinja = nimsave.inj.Ia.*1e-3;
                dat.iinjb = nimsave.inj.Ib.*1e-3;
                dat.iinjc = nimsave.inj.Ic.*1e-3;
                dat.Itor = nimsave.Bprob.surf.Itor_spa_avg.*1e-3;
            catch % psitet
                dat.Itor = psisave.itor;
                dat.ItorTime = psisave.time;
            end
        end
    end
    
    dat.time = out.time;% I DONT KNOW IF THIS IS CORRECT: MAKE SURE IT MATCHES WITH COMPARIRE PLOTS 3
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
    dat.impacts = param.impacts;
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
%    [Phase,Velocity]=VelocityPhase(dat,[.9,1.8],5)
    %cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data');
 %   PhaseVelocity.Phase=Phase;
  %  PhaseVelocity.Velocity=Velocity;
   % save(['Phase' num2str(shots(n)) ], 'PhaseVelocity');
    
    
    cd('T:\IDS\Data Repository'); % \TEMP for changes we dont want to overwrite
    if and(s.spSum, ~s.tSum)
        dat.chSum = s.chSum;
        save(['dat' num2str(shots(n)) 's'], 'dat');
    elseif and(~s.spSum, s.tSum)
        save(['dat' num2str(shots(n)) 't'], 'dat');
    elseif and(s.spSum, s.tSum)
        dat.chSum = s.chSum;
        save(['dat' num2str(shots(n)) 'st'], 'dat');
    else
        save(['dat' int2str(s.sim) num2str(shots(n)) num2str(nBDmodes) ], 'dat') % save 'dat' structure as 'dat<shot><nBDmodes>.mat'
    end
    %end %Try
end
    
profile viewer;
    
    
    
    

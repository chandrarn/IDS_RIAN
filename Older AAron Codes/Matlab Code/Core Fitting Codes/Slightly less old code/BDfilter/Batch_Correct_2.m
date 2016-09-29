% Plot Data, mostly used now to process IDS data for later viewing by
% NIMcompare_1 or similar program.  Modified to incorporate frame adding in
% a spatial or temporal pattern.


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
%140708018,,140709010,140709012:140709017
shots = [129499]; % shot numbers
nBDmodes = [10]; % number of modes to save after BD filtering. Leave blank for no BD filtering
timeBound = []; % time point bounds for BD filtering.  Leave blank to use whole movie
line = [1]; % spectral line number, from longest to shortest wavelength

%Is the file/tree stucture from HIT-SI, or HIT-SI3?
(length(shots(1))>6)
hitsi3=0;
if (length(num2str(shots(1)))>6) 
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
    if ~hitsi3
         HitTree = Tree('hitsi',shots(n));
    else
         HitTree = Tree('hitsi3',shots(n));
    end
    Data = HitTree.getNode('\IDS_IMPACTS');
    NATIVEvalue(Data.getData().data)
    [out, param] = load_shot(shots(n), nBDmodes, timeBound, line, frame_sum, plots, s,hitsi3);
    t = toc;
    disp(['Elapsed time is ' num2str(floor(t/60)) ' minutes, ' num2str(rem(t,60),2) ' seconds']);
    disp(['saving shot ' num2str(shots(n))]);
    
    
    if(~hitsi3)
        %I ASSUME THAT THIS IS HOW TAGS WORK
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
        Data = HitTree.getNode('\IDS_IMPACTS');
        dat.impacts=NATIVEvalue(Data.getData().data);
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
        HitTree = Tree('hitsi3',shots(n));
        Data = HitTree.getNode('\I_INJ_A:RAW');
        iinja=Data.getData().data.getDoubleArray;
        callib = NATIVEvalue(HitTree.getNode('\I_INJ_A:CAL_FACT').getData());
        iinjaTime =  NATIVEvalue(Data.getDimensionAt(0).data());
        dat.iinjaTime = 1e3 * iinjaTime;
        dat.iinja = 1e-3 * iinja * callib;
        Data = HitTree.getNode('\I_INJ_B:RAW');
        iinjb=Data.getData().data.getDoubleArray;
        callib = NATIVEvalue(HitTree.getNode('\I_INJ_B:CAL_FACT').getData());
        iinjbTime =  NATIVEvalue(Data.getDimensionAt(0).data());
        dat.iinjbTime = 1e3 * iinjbTime;
        dat.iinjb = 1e-3 * iinjb * callib;
        Data = HitTree.getNode('\I_INJ_C:RAW');
        iinjc=Data.getData().data.getDoubleArray;
        callib = NATIVEvalue(HitTree.getNode('\I_INJ_C:CAL_FACT').getData());
        iinjcTime =  NATIVEvalue(Data.getDimensionAt(0).data());
        dat.iinjcTime = 1e3 * iinjcTime;
        dat.iinjc = 1e-3 * iinjc * callib;
        Data = HitTree.getNode('\I_TOR_SPAAVG');
        Itor=Data.getData().data.getDoubleArray;
        ItorTime =  NATIVEvalue(Data.getDimensionAt(0).data());
        dat.ItorTime = 1e3 * ItorTime;
        dat.Itor = 1e-3 * Itor;
        Data = HitTree.getNode('\IDS_IMPACTS');
        dat.impacts=NATIVEvalue(Data.getData());
        mdsclose();
    end
    
    dat.time = out.time;%.*1e3;% I DONT KNOW IF THIS IS CORRECT: MAKE SURE IT MATCHES WITH COMPARIRE PLOTS 3
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
%    [Phase,Velocity]=VelocityPhase(dat,[.9,1.8],5)
    %cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data');
 %   PhaseVelocity.Phase=Phase;
  %  PhaseVelocity.Velocity=Velocity;
   % save(['Phase' num2str(shots(n)) ], 'PhaseVelocity');
    
    
    cd('T:\IDS\Data Repository'); % IS THIS RIGHT?
    if and(s.spSum, ~s.tSum)
        dat.chSum = s.chSum;
        save(['dat' num2str(shots(n)) 's'], 'dat');
    elseif and(~s.spSum, s.tSum)
        save(['dat' num2str(shots(n)) 't'], 'dat');
    elseif and(s.spSum, s.tSum)
        dat.chSum = s.chSum;
        save(['dat' num2str(shots(n)) 'st'], 'dat');
    else
        save(['dat' num2str(shots(n)) num2str(nBDmodes)], 'dat') % save 'dat' structure as 'dat<shot><nBDmodes>.mat'
    end
end
    
profile viewer;
    
    
    
    

% Plot Data, mostly used now to process IDS data for later viewing by
% NIMcompare_1 or similar 


%129810,129819,129817,129820
%140708018,
%129810: 14.5
%128585: 53.5
%129213: 68.5
%shots = [151217019:151217026];
%shots = [150624132:program.  Modified to incorporate frame adding in
% a spatial or temporal pattern. Aaron Hossack
% UPDATED TO OBJECT ORIENTED MDSplus BY RIAN CHANDRA, WNTR 2014
% UPDATED TO NEW CINE TYPE BY RIAN CHANDRA, SPR 2014

% NB: Important Input settings Below %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
%close all; 
%clc;
profile on;
import MDSplus.*
Conn=Connection('landau.hit');
global homePath;
homePath = pwd;
homePath=homePath(1:end-21); % grab the base filepath for the git repository
addpath(genpath([homePath '\NewCodes'])); % recursively add files to path
assignin('base','homePath',homePath);


%NOTE: 6/24/15: hardcoded hitsi3 shot to pull from, pdc3 tree is being
%difficult. 


try 
    rmpath('T:\RChandra\A-A-Ron Code\Matlab Code\Core Fitting Codes\');
    addpath('~/IDS/Matlab/') % shot numbers
    end
try
    %addAllThePaths;
end
shots = [160615027:160615028, 160728011:160728018, 160728021:160728024 ,160802008:160802013,160802015:160802020, 160803011:160803012];%,160518015,160518017,160518024,160518029,160518032,160518034:160518036];
nBDmodes = [10]; % number of modes to save after BD filtering. Leave blank for no BD filtering
timeBound = []; % time point bounds for BD filtering.  Leave blank to use whole movie
line = [1:4]; % spectral line number, from longest to shortest wavelength
useTree = 1; % Use ONLY in cases where caliibration was unable to save to the tree
%Is the file/tree stucture from HIT-SI, or HIT-SI3?
% (length(shots(1))>6)

if (length(num2str(shots(1))) > 8) 
    hitsi3 = 1;
else
    hitsi3 = 0;
end
    
label1 = 'Impact Parameter [cm]'; % Toroidal Section
% label1 = 'Major Radius [cm]'; % Poloidal Section, fibers facing each other
label2 = 'Major Radius [cm]'; % Poloidal Section
% label1 = 'Displacement [cm]'; % Orthogonal Section

frame_sum = 1; % number of frames to sum (typically 1)
plots = 0; % makes plots for debugging

linuxDataPath = '/home/aaron/IDS/IDSdata/';
windowsDataPath = 'T:\IDS\Data Repository\';

s.sim = [];  % tag indicating type of simulation
              % 2 = NIMROD --- weighted by n^2 ---  NO temperature  --- 7x timeSum
              % 4 = NIMROD --- weighted by n^2 ---  NO temperature  --- no timeSum
              % 5 = PSI-TET --- weighted by n^2 --- with temperature --- 7x timeSum
              % 6 = NIMROD --- weighted by n^2 --- with temperature --- 7x timeSum
              % 7 = PSI-TET -- weighted by n^2 --- with temperature --- no timeSum
              % 8 = NIMROD --- weighted by n^2 --- with temperature --- no timeSum
              % LEAVE EMPTY [] IF NOT DOING SIMULATION

s.spSum = 0; % add chords together for better statistics
%     s.chSum = [10:14; 23:27]; % channels to add together.  The columns of each row are added together.
s.chSum = [9:11; 12:14; 16:18; 19:21; 23:25; 26:28];

s.tSum = 0; % add frames together in a temporal pattern
s.sumPatt = [50, 12, 200]; % [time point to start from, length of pattern, end]

%addpath('T:\RChandra\A-A-Ron Code\Matlab Code\Core Fitting Codes');
%addpath('T:\RChandra\A-A-Ron Code\General Matlab');
%addpath('T:\IDS\Data Analysis');
for n = 1:length(shots)
    tic
%     if hitsi3
%         HitTree = Tree('analysis3',shots(n));
%     else
%         HitTree = Tree('analysis',shots(n));
%     end
%     Data = Conn.get('\IDS_IMPACTS');
%     NATIVEvalue(Data.getData().data);
%     cd('T:\IDS\Data Analysis');
    [out, param] = load_shot(shots(n), nBDmodes, timeBound, line, frame_sum, plots, s, hitsi3, useTree);
    t = toc;
    disp(['Elapsed time is ' num2str(floor(t/60)) ' minutes, ' num2str(rem(t,60),2) ' seconds']);
    if ~isempty(s.sim)        
        disp(['saving shot ' num2str(s.sim) num2str(shots(n))]);
    else
        disp(['saving shot ' num2str(shots(n))]);
    end
    
    if(~hitsi3)
        %I ASSUME THAT THIS IS HOW TAGS WORK
        HitTree = Tree('hitsi',shots(n));
        Data = Conn.get('\I_INJ_X');
        iinjx=Data.getData().data.getDoubleArray;
        %[iinjxTime, dummy, iinjx] = Data.getData();
        iinjxTime =  NATIVEvalue(Data.getDimensionAt(0).data());
        dat(1).iinjxTime = 1e3 * iinjxTime;
        dat(1).iinjx = 1e-3 * iinjx;
        Data = Conn.get('\I_INJ_Y');
        iinjy = Data.getData().data.getDoubleArray;
        iinjyTime = NATIVEvalue(Data.getDimensionAt(0).data());
        %[iinjyTime, dummy, iinjy] = Data.getData();
        dat(1).iinjyTime = 1e3 * iinjyTime;
        dat(1).iinjy = 1e-3 * iinjy;
        Data = Conn.get('\I_TOR_SPAAVG');
        %[ItorTime, dummy, Itor]=Data.getData();
        Itor = Data.getData().data.getDoubleArray;
        ItorTime =  NATIVEvalue(Data.getDimensionAt(0).data());
        dat(1).ItorTime = 1e3 * ItorTime;
        dat(1).Itor = 1e-3 * Itor;
       
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
        Conn.openTree('hitsi3',shots(n));%Tree('hitsi3',shots(n));
        dat(1).iinja = 1e-3 * (Conn.get('\I_INJ_A').getFloatArray);
        dat(1).iinjaTime = 1e3 *  Conn.get('dim_of(\I_INJ_A)').getFloatArray;
        dat(1).iinjb = 1e-3 *Conn.get('\I_INJ_B').getFloatArray;
        dat(1).iinjbTime = 1e3 *  Conn.get('dim_of(\I_INJ_B)').getFloatArray;
        dat(1).iinjc = 1e-3 *Conn.get('\I_INJ_C').getFloatArray;
        dat(1).iinjcTime = 1e3 *  Conn.get('dim_of(\I_INJ_C)').getFloatArray;
        Itor = Conn.get('\I_TOR_SPAAVG').getFloatArray;
        ItorTime =  Conn.get('dim_of(\I_TOR_SPAAVG)').getFloatArray;
        dat(1).ItorTime = 1e3 * ItorTime;
        dat(1).Itor = 1e-3 * Itor;
        
        mdsclose();
    end
    
    for nn = 1:length(line)
        dat(nn).temp = out(nn).temp;
        dat(nn).vel = 1e-3 * out(nn).vel; % save in km/s
        dat(nn).int = out(nn).int;
        dat(nn).fit_par = out(nn).fit_par; % Curve fit parameters [time, channel, 1:6]
        dat(nn).bounds = out(nn).bounds; % Bounds for curve fitting grid, (time) x (channel) x (x1 x2 y1 y2)
        dat(nn).guesses = out(nn).guesses; % initial guess parameters, same size as fit_par
        if param.calcError
            dat(nn).dparam = out(nn).dparams; % parameter uncertainty structure
            dat(nn).stddev = out(nn).stddev; % standard deviation of every data point
            dat(nn).intU = out(nn).intU; % upper error bar for intensity
            dat(nn).intL = out(nn).intL;
            dat(nn).tempU = out(nn).tempU; % upper error bar for temperature
            dat(nn).tempL = out(nn).tempL;
            dat(nn).velU = 1e-3 * out(nn).velU; % upper error bar for velocity
            dat(nn).velL = 1e-3 * out(nn).velL;
        end
        if exist('out.resnorm', 'var') % DNE if using LM curve fitting method
            dat(nn).residual = sqrt(out(nn).resnorm); % residual
        end
        if exist('out.exp', 'var')
            dat(nn).exp = out(nn).exp;
            dat(nn).corrected = out(nn).corrected;
        end
    end
    
    dat(1).time = out(1).time; 
    dat(1).raw = out(1).data; % raw data points
    dat(1).peaks = param.peaks(:, 1);
    dat(1).param = param; % param structure
    dat(1).shotRef = param.shotRef; % real shot number
    dat(1).impacts = param.impacts;
    dat(1).title = ['Shot ' num2str(dat(1).shotRef)];
    dat(1).label1 = label1;
    dat(1).label2 = label2;
  
    %Average chord velocity by injector phase. WHERE TO SAVE?
%    [Phase,Velocity]=VelocityPhase(dat,[.9,1.8],5)
    %cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data');
 %   PhaseVelocity.Phase=Phase;
  %  PhaseVelocity.Velocity=Velocity;
   % save(['Phase' num2str(shots(n)) ], 'PhaseVelocity');
    
    
%     cd('T:\IDS\Data Repository\TEMP'); % Temp for changes we dont want to overwrite
%     if and(s.spSum, ~s.tSum)
%         dat.chSum = s.chSum;
%         save(['dat' num2str(shots(n)) 's'], 'dat');
%     elseif and(~s.spSum, s.tSum)
%         save(['dat' num2str(shots(n)) 't'], 'dat');
%     elseif and(s.spSum, s.tSum)
%         dat.chSum = s.chSum;
%         save(['dat' num2str(shots(n)) 'st'], 'dat');
%     else
try
    save([linuxDataPath 'dat' num2str(s.sim) num2str(shots(n)) num2str(nBDmodes)], 'dat') % save 'dat' structure as 'dat<shot><nBDmodes>.mat'
catch
    save([windowsDataPath 'dat' num2str(s.sim) num2str(shots(n)) num2str(nBDmodes)], 'dat') % save 'dat' structure as 'dat<shot><nBDmodes>.mat'
end
%     end
end
    
profile viewer;
    
    
    
    

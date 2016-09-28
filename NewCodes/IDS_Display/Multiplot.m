% Multiplot
% Note: need to manually save the figure as .eps, otherwise it saves in
% black and white. You know, because fuck all.
%close all; 
clear all;% clc;
addpath('~/IDS/Matlab/');
addpath('T:\RChandra\Sine_fit\');
addpath('T:\IDS\General Matlab\');
%addAllThePaths;
lines = {'O II', 'C III', 'O II','C III'};

%% Input Settings

%% IDS DATA
% in(1).shot = 160601022;150625998;
% in(1).line =3; % line # NB: 1 is C III, 2 is O II, 3 is C III !
% in(1).legend = [num2str(in(1).shot) ' ' lines{in(1).line}];
% in(1).color = {'r';[225,105,0]./255};
% in(1).style = '-';
% in(1).error = 0; % 1 / 0 for errorbars
% in(1).velShift = -5; % SHIFT VELOCITY
% in(1).intScale = 1; % scale factor for intensity
% in(1).timeShift = 0; % ms, shift time base
% in(1).timeScale = 1e-3; % scale timebase to put into ms
% in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
% in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(1).doubleplot = [1];[1:23; 24,26:47]; % plot coorespoinding impacts
% in(1).fftPlot = [1]; % FFT of signal, n frequencies
% 129499 
% note: need to change .*1e-6 to .*1e-3 in sinefit
in(1).shot = 129499;%150625998;
in(1).line = 2; % line # NB: 1 is C III, 2 is O II, 3 is C III !
in(1).legend = [num2str(in(1).shot) ' O II'];
in(1).color = {'r';[225,105,0]./255};
in(1).style = '-';
in(1).error = 0; % 1 / 0 for errorbars
in(1).velShift = 5; % SHIFT VELOCITY
in(1).intScale = 1; % scale factor for intensity
in(1).timeShift = 0; % ms, shift time base
in(1).timeScale = 1;1e-3; % scale timebase to put into ms
in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
in(1).doubleplot = [];[1:23; 24,26:47]; % plot coorespoinding impacts
in(1).fftPlot = [1]; % FFT of signal, n frequencies
% 
% in(2).shot = 151217026;%150625998;
% in(2).line = 2; % line # NB: 1 is O II, 2 is C III !
% in(2).legend = [num2str(in(1).shot) ' O II'];
% in(2).color = {'b';[0,109,200]./255};
% in(2).style = '-';
% in(2).error = 0; % 1 / 0 for errorbars
% in(2).velShift = 5; % SHIFT VELOCITY
% in(2).intScale = 1; % scale factor for intensity
% in(2).timeShift = 0; % ms, shift time base
% in(2).timeScale = 1e-3; % scale timebase to put into ms
% in(2).injTimeScale = 1e0;1e-3; % scale the injector time to ms
% in(2).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(2).doubleplot = [1:23; 24,26:47]; % plot coorespoinding impacts
% in(2).fftPlot = [1]; % FFT of signal, n frequencies

% in(3).shot = 151217026;%150625998;
% in(3).line = 3; % line # NB: 1 is O II, 2 is C III !
% in(3).legend = [num2str(in(1).shot) ' O II'];
% in(3).color = {[0,128,102]./255;[117,186,102]./255};
% in(3).style = '-';
% in(3).error = 0; % 1 / 0 for errorbars
% in(3).velShift = 5; % SHIFT VELOCITY
% in(3).intScale = 1; % scale factor for intensity
% in(3).timeShift = 0; % ms, shift time base
% in(3).timeScale = 1e-3; % scale timebase to put into ms
% in(3).injTimeScale = 1e0;1e-3; % scale the injector time to ms
% in(3).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(3).doubleplot = [1:23; 24,26:47]; % plot coorespoinding impacts
% in(3).fftPlot = [1]; % FFT of signal, n frequencies
% 
% in(2).shot = 12949610;
% in(2).line = 2; % line #
% in(2).legend = '129496 C III';
% in(2).color = 'r';
% in(2).style = '-';
% in(2).error = 1; % 1 / 0 for errorbars
% in(2).velShift = 0; % SHIFT VELOCITY
% in(2).intScale = 1; % scale factor for intensity
% in(2).timeShift = 0; % ms, shift time base

%% 129530, NIMROD, PSI-TET
% in(1).shot = 12953010;
% in(1).line = 2; % line # NB: 1 is O II, 2 is C III !
% in(1).legend = '129530 C III';
% in(1).color = 'b';
% in(1).style = '-';
% in(1).error = 1; % 1 / 0 for errorbars
% in(1).velShift = 0; % SHIFT VELOCITY
% in(1).intScale = 1; % scale factor for intensity
% in(1).timeShift = 0; % ms, shift time base

%% 129499, NIMROD, PSI-TET
% in(1).shot = 12949910;
% in(1).line = 2; % line # NB: 1 is O II, 2 is C III !
% in(1).legend = '129499 C III';
% in(1).color = 'b';
% in(1).style = '-';
% in(1).error = 1; % 1 / 0 for errorbars
% in(1).velShift = 0; % SHIFT VELOCITY
% in(1).intScale = 1; % scale factor for intensity
% in(1).timeShift = 0; % ms, shift time base
% 
% in(2).shot = 6129499;
% in(2).line = 1; % line #
% in(2).legend = 'NIMROD';
% in(2).color = 'r';
% in(2).style = '-';
% in(2).error = 0; % 1 / 0 for errorbars
% in(2).velShift = 0; % SHIFT VELOCITY
% in(2).intScale = 0.35; % scale factor for intensity
% in(2).timeShift = 0.9377; % ms, shift time base
% 
% in(3).shot = 5129499;
% in(3).line = 1; % line #
% in(3).legend = 'PSI-TET';
% in(3).color = [12/255 117/255 0];
% in(3).style = '-';
% in(3).error = 0; % 1 / 0 for errorbars
% in(3).velShift = 0; % SHIFT VELOCITY
% in(3).intScale = 0.35; % scale factor for intensity
% in(3).timeShift = 1.335; % ms, shift time base

%% NIMROD temperature details
% in(1).shot = 6129499;
% in(1).line = 1; % line #
% in(1).legend = 'w/ Temp.';
% in(1).color = 'r';
% in(1).style = '-';
% in(1).error = 0; % 1 / 0 for errorbars
% in(1).velShift = 0; % SHIFT VELOCITY
% in(1).intScale = 0.35; % scale factor for intensity
% in(1).timeShift = 0.9377; % ms, shift time base
% 
% in(2).shot = 2129499;
% in(2).line = 1; % line #
% in(2).legend = 'NO Temp.';
% in(2).color = 'm';
% in(2).style = '-';
% in(2).error = 0; % 1 / 0 for errorbars
% in(2).velShift = 0; % SHIFT VELOCITY
% in(2).intScale = 0.35; % scale factor for intensity
% in(2).timeShift = 0.9377; % ms, shift time base
% 
% in(3).shot = 8129499;
% in(3).line = 1; % line #
% in(3).legend = 'NO Time Avg.';
% in(3).color = [0.5 0 0];
% in(3).style = '--';
% in(3).error = 0; % 1 / 0 for errorbars
% in(3).velShift = 0; % SHIFT VELOCITY
% in(3).intScale = 0.35; % scale factor for intensity
% in(3).timeShift = 0.9377; % ms, shift time base
% 
% in(4).shot = 62129499;
% in(4).line = 1; % line #
% in(4).legend = 'Vel. Subtracted';
% in(4).color = [0.6 0 0.4];
% in(4).style = '-.';
% in(4).error = 0; % 1 / 0 for errorbars
% in(4).velShift = 0; % SHIFT VELOCITY
% in(4).intScale = 0.35; % scale factor for intensity
% in(4).timeShift = 0.9377; % ms, shift time base

%% PSI-TET temperature details
% in(1).shot = 5129499;
% in(1).line = 1; % line #
% in(1).legend = 'w/ Temp.';
% in(1).color = [12/255 117/255 0];
% in(1).style = '-';
% in(1).error = 0; % 1 / 0 for errorbars
% in(1).velShift = 0; % SHIFT VELOCITY
% in(1).intScale = 0.35; % scale factor for intensity
% in(1).timeShift = 1.335; % ms, shift time base
% 
% in(2).shot = 3129499;
% in(2).line = 1; % line #
% in(2).legend = 'NO Temp.';
% in(2).color = 'g';
% in(2).style = '-';
% in(2).error = 0; % 1 / 0 for errorbars
% in(2).velShift = 0; % SHIFT VELOCITY
% in(2).intScale = 0.35; % scale factor for intensity
% in(2).timeShift = 1.335; % ms, shift time base
% 
% in(3).shot = 7129499;
% in(3).line = 1; % line #
% in(3).legend = 'NO Time Avg.';
% in(3).color = [6/255 59/255 0];
% in(3).style = '--';
% in(3).error = 0; % 1 / 0 for errorbars
% in(3).velShift = 0; % SHIFT VELOCITY
% in(3).intScale = 0.35; % scale factor for intensity
% in(3).timeShift = 1.335; % ms, shift time base
% 
% in(4).shot = 53129499;
% in(4).line = 1; % line #
% in(4).legend = 'Vel. Broadening Subtracted';
% in(4).color = [6/255 59/255 0];
% in(4).style = '-.';
% in(4).error = 0; % 1 / 0 for errorbars
% in(4).velShift = 0; % SHIFT VELOCITY
% in(4).intScale = 0.35; % scale factor for intensity
% in(4).timeShift = 1.335; % ms, shift time base

%%
saving = 0;
plotCurrents = 1;
plotAverages = 0;
compactCurrents = 1;
plotTor = 1; % The line plot will show the differenve between fibers

if isempty(in(1).fftPlot)
    Analysis=1;
else
    Analysis = 2; % Analyze torroidal flow, Amplitude, phasing. Replaces "Averages"
end

 plotType = 1; % Velocity
% plotType = 2; % Temperature
% plotType = 3; % Intensity

timebound = [1.1, 2.0]; % [ms]
% timebound = [1.2 2.3]; % [ms]

CutPow  = .4; % FFT Power Cuttoff Percentage. 

saveFile = ['T:\IDS\Analysis Repository\' num2str(in(1).shot)];

% Set these: ( add :2: if you only want every other line
% chan_range = 50:58;
% defaults
chan_ranget = [7:26];
chan_rangep = [43:62];
% 151217026 values:
if in(1).shot == 151217026
    chan_ranget = [13:26];[3:62];[8:26]; % toroidal, mohawk port in midplane
    chan_rangep = [49:62];[44:62]; % poloidal
    timebound=[.8,2.0];
elseif in(1).shot == 151217024
    chan_ranget = [11:26];
    chan_rangep = [47:62];
    timebound=[.8,2.0];
elseif in(1).shot == 151217025
    chan_ranget = [13:26];
    chan_rangep = [49:62];
    timebound=[.8,2.0];
elseif in(1).shot == 151217020 || in(1).shot == 151217021 || in(1).shot == 151217019 || in(1).shot == 151217022|| in(1).shot == 151217016
    chan_ranget = [15:25];
    chan_rangep = [51:61];
elseif in(1).shot == 129499
    chan_ranget = [4:27];
    chan_rangep = []; % we dont want the axial port
    timebound = [1.4 2.0];
elseif in(1).shot == 160525016
    timebound = [1.5,2.0];
elseif in(1).shot == 160525017
    timebound = [1.35,2.0];
elseif in(1).shot >= 160526030 && in(1).shot <= 160526038
    timebound = [1.35,2.0];
    chan_ranget = [12:23];
    chan_rangep = [48:59];
end
chan_range = [chan_ranget, chan_rangep];
xlim = [0,50];
%chan_range = chan_rangep
%chan_range = [1:68];
%chan_range = [8:2:24];
% chan_range = [43:53];

velSpace = 15; % km/s
intSpace = 20; % arb.u.
tempSpace = 20; % eV

%% set up figure
fntsz = 14;
lnwdth = 1.5;
errWdth = 500; % errorbar width setting

S = get(0, 'ScreenSize');
analysisHeight = S(4) - 1;
figureWidth = (S(3) - 12)/2.25;
colors = [ 0 0 1; 12/255 117/255 0; 1 0 0 ];
h = figure('Visible', 'on', 'Name', ['MULTIPLOT-Lines: ' num2str(in(1).line)], 'Position',...
    [5, 35, figureWidth, analysisHeight], 'Color', [1 1 1]);

if ~compactCurrents
    ax = axes('Parent', h, 'Position', [0.075, 0.08, 0.8, 0.85]);
else
    ax = axes('Parent', h, 'Position', [0.075, 0.28, 0.8, 0.65]);
end

hold on;

if plotAverages || Analysis ==1
    h2 = figure('Visible', 'on', 'Name', ['MULTIPLOT-Toroidal Flow: ' num2str(in(1).line)], 'Position',...
        [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
    ax2 = axes('Parent', h2, 'Position', [0.075, 0.15, 0.85, 0.75]);
    hold on;
    grid on;
end
if Analysis==1 && ~isempty(in(1).fftPlot)
     h3 = figure('Visible', 'on', 'Name', ['MULTIPLOT-Phase: ' num2str(in(1).line)], 'Position',...
        [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
    ax3 = axes('Parent', h3, 'Position', [0.075, 0.15, 0.85, 0.75]); hold on; grid on;
    h4 = figure('Visible', 'on', 'Name', ['MULTIPLOT-Displacement: ' num2str(in(1).line)], 'Position',...
        [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
    ax4 = axes('Parent', h4, 'Position', [0.075, 0.15, 0.85, 0.75]); hold on; grid on;
     h5 = figure('Visible', 'on', 'Name', ['MULTIPLOT-FlowRanges: ' num2str(in(1).line)], 'Position',...
        [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
    ax5 = axes('Parent', h5, 'Position', [0.075, 0.15, 0.85, 0.75]);
    hold on;
    grid on;
elseif Analysis==2 && ~isempty(in(1).fftPlot)
    h2=figure;
    h6 = figure('Visible', 'on', 'Name', ['MULTIPLOT-Analysis: ' in(1).legend], 'Position',...
    [5, 1, figureWidth-300, analysisHeight], 'Color', [1 1 1]);
    ax6 = axes('Parent',h6,'Position',[0.15,.08,.8,.25]); hold on; grid on;box on;
    ax7 = axes('Parent',h6,'Position',[0.15,.38,.8,.25]); hold on; grid on;box on;
    ax8 = axes('Parent',h6,'Position',[0.15,.68,.8,.25]); hold on; grid on;box on;
    txt=text(10,550,['Analysis: ' in(1).legend],'fontweight','bold','fontsize',13);
    title(ax8,'Phase');
    title(ax7,'Toroidal Flow');
    title(ax6,'Average Displacement');
    h7 = figure('Visible', 'on', 'Name', ['MULTIPLOT-FFT Spectrum: ' num2str(in(1).line)], 'Position',...
        [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
    ax9 = axes('Parent', h7, 'Position', [0.075, 0.15, 0.85, 0.75]); hold on; grid on;box on;
%     h3=figure;
%     h4=figure;
%     h5=figure;
end
figure(h); % make first figure current

%% Load and Plot Data
clear param
for n = 1:length(in)
    addpath('T:\IDS\Data Repository');
    
    if ~(exist('dat','var') && strcmp(dat(1).title,['Shot ' num2str(in(n).shot)])) % speed up runtime if data is already loaded
        load(['dat' num2str(in(n).shot) '10.mat']); % Real HIT-SI Data
    end

    dat = trimRange(dat, chan_range, 0,timebound.*(1/in(n).timeScale),[]);
    Itor(:,n) = dat(1).Itor;
    
    switch plotType
        case 1
            if ~isempty(in(n).doubleplot) && isempty(in(1).fftPlot)
                %data(1:length(dat(1).time),:) = dat(in(n).line).vel(:,in(n).doubleplot(1,:));
                %data(length(dat(1).time)+1:2*length(dat(1).time),:) = ...
                %    dat(in(n).line).vel(:,in(n).doubleplot(2,:));
                data(1:length(dat(1).time),:) = dat(in(n).line).vel(:,1:(length(dat(1).impacts))/2);
                data(length(dat(1).time)+1:2*length(dat(1).time),:) = ...
                   dat(in(n).line).vel(:,(length(dat(1).impacts)/2)+1:end);
            elseif ~isempty(in(1).fftPlot)
                dat(in(n).line).vel = averageNans(dat(in(n).line).vel); % remove nans
                
               
                if ~isempty(in(n).doubleplot)
                    % Find where each fiber bundle begins and ends.
                    doubleplot(1,:) = 1:(length(dat(1).impacts))/2;
                    doubleplot(2,:) = (length(dat(1).impacts)/2)+1:length(dat(1).impacts);

                    % initialize Sine_Fit parameters
                    param(:,1,n) = dat(1).impacts(doubleplot(1,:));
                    param(:,6,n) = dat(1).impacts(doubleplot(2,:));
                else
                    doubleplot(1,:) = 1:length(dat(1).impacts);
                    param(:,1,n) = dat(1).impacts(doubleplot(1,:));
                end
                 % Param: impacts offset amplitude phase, frequency
                 if ~exist('param','var')
                        param = zeros(length(doubleplot),10,length(in));
                 end
                display('Computing FFT');
                pRel = zeros(length(doubleplot),2);
                for i = 1:length(doubleplot)
                    % using fftf method
%                    offset = mean(dat(in(n).line).vel(:,in(n).doubleplot(1,i)));
%                    [data(1:length(dat(1).time),i), f, y, y2] = ...
%                        fftf(dat(1).time.*1e-6, dat(in(n).line).vel( ...
%                        :,in(n).doubleplot(1,i)) - offset, 50e3,3,10e3,1);
%                    data(1:length(dat(1).time),i) = ...
%                    data(1:length(dat(1).time),i).* abs(max(y2)); % make real
%                    data(1:length(dat(1).time),i) = ...
%                        data(1:length(dat(1).time),i) + offset; % replace offset
%                    pause(1);
%                    offset = mean(dat(in(n).line).vel(:,in(n).doubleplot(2,i)));
%                    [data(length(dat(1).time)+1:2*length(dat(1).time),i), f, y, y2] = ...
%                        fftf(dat(1).time.*1e-6, dat(in(n).line).vel( ...
%                        :,in(n).doubleplot(2,i))- offset, 50e3,3,10e3,1);
%                    data(length(dat(1).time)+1:2*length(dat(1).time),i) = ...
%                        data(length(dat(1).time)+1:2*length(dat(1).time),i).* abs(max(y2));
%                    data(length(dat(1).time)+1:2*length(dat(1).time),i) = ...
%                        data(length(dat(1).time)+1:2*length(dat(1).time),i) + offset;
%                    pause(1);
                    % using SineFit method
                    %try
                        signal = dat(in(n).line).vel(:,doubleplot(1,i))+in(n).velShift;
                        Fsamp = 1/(mean(diff(dat(1).time.*(in(n).timeScale.*1e-3))));
                        offset = nanmean(signal);
                        amp = max(signal)-offset;
                        freq = 14500;
                        phase = pi/2;
                        [param(i,2:5,n),data(1:length(dat(1).time),i)] = sine_fit( ...
                            dat(1).time'.*(in(n).timeScale.*1e-3),signal',[nan,nan,nan,freq], ...
                            [offset,amp,phase,freq],0);
                        if param(i,3,n)<0 % 180degree phase
                            param(i,3,n)=-param(i,3,n);
                            param(i,4,n)=param(i,4,n)+pi;
                            disp([' WARNING: NEGATIVE AMPLITUDE @ n=' num2str(n) ', i=' num2str(i)]);
                        end
                        
                        
                        % Calculate fft reconstruction validity
                        harm1 = bandpower(signal-offset,Fsamp,[13,17].*1e3);
                        harm2 = bandpower(signal-offset,Fsamp,[26,31].*1e3);
                        ptot =  bandpower(signal-offset,Fsamp,[0,Fsamp*(length(signal)-2)/(2*length(signal)) ]); % Nyquist
                        pRel(i,1) = (harm1+harm2)/ptot;
                        % if the fit isnt valid, dont plot it
                        try data(1:length(dat(1).time),i.*(pRel(i,1)<CutPow))=signal;end
                        
                        %pause(1);
                        if in(n).doubleplot
                        signal = dat(in(n).line).vel(:,doubleplot(2,i))+in(n).velShift;
                        offset = mean(signal);
                        amp = max(signal)-offset;
                        [param(i,7:10,n),data(length(dat(1).time)+1:2*length(dat(1).time),i)] = ...
                            sine_fit(dat(1).time'.*(in(n).timeScale.*1e-3),signal',[nan,nan,nan,freq], ...
                            [offset,amp,phase,freq],0);
                        %pause(1);
                        if param(i,8,n)<0 % 180degree phase
                            param(i,8,n)=-param(i,8,n);
                            param(i,9,n)=param(i,9,n)+pi;
                            disp([' WARNING: NEGATIVE AMPLITUDE @ n=' num2str(n) ', i=' num2str(i)]);
                        end
                        % Calculate fft reconstruction validity
                        harm1 = bandpower(signal-offset,Fsamp,[13,17].*1e3);
                        harm2 = bandpower(signal-offset,Fsamp,[26,31].*1e3);
                        ptot =  bandpower(signal-offset,Fsamp,[0,Fsamp*(length(signal)-2)/(2*length(signal))]);
                        pRel(i,2) = (harm1+harm2)/ptot;
                        % if the fit isnt valid, dont plot it
                        try data(length(dat(1).time)+1:2*length(dat(1).time)...
                                ,i.*(pRel(i,2)<CutPow))=signal;end
                        end
                    
                    %catch
                    %    display(['halted at ' num2str(i)]);
                    %end
                end
            else
                data = dat(in(n).line).vel + in(n).velShift;
            end
            titles = 'Velocities';
            offset = velSpace;
            units = 'km/s';
            sidebar = [num2str(offset) ' km/s per division'];
            if in(n).error
                errorL = dat(in(n).line).velL;
                errorU = dat(in(n).line).velU;
            end
        case 2
            % Check Double Plotting
             if ~isempty(in(n).doubleplot) && isempty(in(1).fftPlot)
                doubleplot(1,:) = 1:(length(dat(1).impacts))/2;
                doubleplot(2,:) = (length(dat(1).impacts)/2)+1:length(dat(1).impacts);
                data(1:length(dat(1).time),:) = dat(in(n).line).temp(:,doubleplot(1,:));
                data(length(dat(1).time)+1:2*length(dat(1).time),:) = ...
                    dat(in(n).line).temp(:,doubleplot(2,:));
             else
            data = dat(in(n).line).temp;
             end
             
            titles = 'Temperatures';
            offset = tempSpace;
            units = 'eV';
            sidebar = [num2str(offset) ' eV per division'];
            if in(n).error
                errorL = dat(in(n).line).tempL;
                errorU = dat(in(n).line).tempU;
            end
        case 3
            data = in(n).intScale * dat(in(n).line).int;
            titles = 'Intensities';
            sidebar = ['Arb.'];
            offset = intSpace;
            units = 'Arb.';
            if in(n).error
                errorL = dat(in(n).line).intL;
                errorU = dat(in(n).line).intU;
            end
    end
    
    
    
    if in(n).doubleplot
         time(1:length(dat(1).time)) = dat(1).time.*in(1).timeScale + in(n).timeShift;
         time(length(dat(1).time)+1:2*length(dat(1).time)) =  time(1:length(dat(1).time));
%              dat(1).time.*in(1).timeScale + in(n).timeShift;
    else
        time = dat(1).time.*in(1).timeScale + in(n).timeShift;
    end

    %% Averages
    if plotAverages
        figure(h2) % make current
        
        % find index corresponding to time bounds
        nTimeLim(1) = find(dat(1).time >= timebound(1), 1);
        nTimeLim(2) = find(dat(1).time <= timebound(end), 1, 'last');
%         nTimeLim(2) = nTimeLim(2) - 1; % the above command goes one too far
        
        % Calculate all data
        for m = 1:size(data, 2);
            selection = data(nTimeLim(1):nTimeLim(2), m);
            dataAvg(m) = mean(selection(~isnan(selection)));
            dataStd(m) = std(selection(~isnan(selection)));
        end
        
        % Plot Data
        t2(n) = errorbar(dat(1).impacts, dataAvg, dataStd, 'color', in(n).color, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
        
        figure(h) % make other current
    end
    
    %% Analysis
    if Analysis
        if plotType ==1
            figure(h2)
            lines=1;
            if ~isempty(in(n).doubleplot)
                lines=2;
            end


            for i = 1:lines
                if in(n).fftPlot % if computing the fft, use offset, amp 
                    % Toroidal Current
                    dataAvg(:,i) = param(:,2+5*(i-1),n);
                    dataStd(:,i) = param(:,3+5*(i-1),n);
                    
                    % Phase
                    %dataPhase(:,i) = param(:,4+5*(i-1),n);
                    % test automatic phase finding
                    for j = 1: size(param,1) % for all impacts
                        dataPhase(j,i) = find_Phase(param(j,1+(5*(i-1)): 5 + 5*(i-1)));
                    end
                    
                    %dataPhase(:,i) = mod(dataPhase(:,i),2*pi); % 2Pi Periodicity
                    
                    %% HARDCODED DATA
                    if in(n).line==1 && strcmp(dat(1).title,'Shot 151217026'); % FOR 151217026 LINE 1
                       % dataPhase(11-1,1) = 6.14; % real weak signal, this is necessary because reasons
                        %dataPhase(21,1)= 3.5; % When signal is this weak, make it match the other fiber
                    elseif in(n).line==3 && strcmp(dat(1).title,'Shot 151217026')
                        if i==1
                            %dataPhase(:,1) = dataPhase(:,1)+2*pi;
                          %  dataPhase(2-1,1) = 2*pi;
                          %  dataPhase(4-1,1) = 5.8;
                        end
                    elseif in(n).line==2 && strcmp(dat(1).title,'Shot 151217026')
                        if i==2
                         %   dataPhase(7-1,2)=5;
                        end
                    elseif in(n).line==1 && strcmp(dat(1).title,'Shot 151217021')
                       % dataPhase(2,1)=0;
                        dataPhase(1,1)=3.2;
                        
                    end
                    %dataPhase(:,i) = mod(dataPhase(:,i),5.43); % Minimum resolution, 151217026
                    
                    % shift lines if a periodicity jump occurs
                    for j = 2:size(data, 2)
                        if (dataPhase(j,i) - dataPhase(j-1,i))>(pi)
                            dataPhase(j,i) = dataPhase(j,i)-(2*pi);
                            disp(['JUMP: -2Pi, Line: ' num2str(i) ' Impact: ' num2str(dat(1).impacts(j))]);
                        elseif (dataPhase(j,i) - dataPhase(j-1,i))<(-pi)
                           disp(['JUMP: +2Pi, Line: ' num2str(i) ' Impact: ' num2str(dat(1).impacts(j))]);
                            dataPhase(j,i) = dataPhase(j,i)+(2*pi);
                        end
                    end
                    
                    
                    % half of integral of one half period is radius of motion
                    dataDispl(:,i) = 2*param(:,3+5*(i-1),n).*(1./param(:,5+5*(i-1),n))./(2*pi)./2 .*1e5;

                else % if no FFT
                    i
                     % find index corresponding to time bounds
                    nTimeLim(1) = find(dat(1).time.*in(n).timeScale >= timebound(1), 1);
                    nTimeLim(2) = find(dat(1).time.*in(n).timeScale <= timebound(end), 1, 'last');
            %         nTimeLim(2) = nTimeLim(2) - 1; % the above command goes one too far

                    % Calculate all data
                    for m = 1:size(data, 2);
                        selection = data((nTimeLim(1):nTimeLim(2))+(i-1)*size(data,1)/2, m);
                        dataAvg(m,i) = mean(selection(~isnan(selection)));
                        dataStd(m,i) = std(selection(~isnan(selection)));
                        dataDispl(m,i) = 2*dataStd(m,i).*(1/14500)./(2*pi)./2 .*1e5;

                    end
                end

                 % Plot Data
                 hold on;
                if Analysis==2
                    %for j=1:length(dataAvg(:,i));try dataAvg(j.*(pRel(j,i)<.5),i)=NaN;end;end
                    t2(n) = errorbar(ax7,dat(1).impacts(1:size(data,2)), dataAvg(:,i), dataStd(:,i), 'color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
                elseif Analysis==1
                    t2(n) = errorbar(dat(1).impacts(1:size(data,2)), dataAvg(:,i), dataStd(:,i), 'color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
                end
                if in(n).fftPlot
                    
                    if Analysis==1
                        figure(h4);
                        hold on;
                        plot(dat(1).impacts(1:size(data,2)),dataDispl(:,i),'color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
                        ylabel('Average Radius of Displacement [cm]');set(gca,'ylim',[0,8]);
                        figure(h2);
                    elseif Analysis==2
                        for j=1:length(dataDispl(:,i));try dataDispl(j.*(pRel(j,i)<CutPow),i)=NaN;end;end
                         plot(ax6,dat(1).impacts(1:size(data,2)),dataDispl(:,i),'-*','color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
                        ylabel(ax6,'[cm]');set(ax6,'ylim',[0,8]);
                        xlabel(ax6,'Impacts [cm]'); set(ax6,'xlim',xlim);
                    end
                end
            end
            if in(n).fftPlot & in(n).doubleplot==1
                for i = 1:size(data, 2)
                    cycle1 = (param(i,3,n)*sin(param(i,4,n)+(2*pi).*(0:(1/100):1))+param(i,2,n));
                    cycle2 = (param(i,8,n)*sin(param(i,9,n)+(2*pi).*(0:(1/100):1))+param(i,7,n));
                    maxDispl(i) = max(cycle1-cycle2);
                    minDispl(i) = min(cycle1-cycle2);
                end
                
                L=-(dataAvg(:,1)-dataAvg(:,2))+minDispl';
                U=-(dataAvg(:,1)-dataAvg(:,2))-maxDispl';
                if Analysis==1
                    figure(h5)
                    errorbar(dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2)),L,U,'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style);
                    set(gca,'ylim',[-2,20]);
                    
                    figure(h3);
                elseif Analysis==2
                    %for k=1:2;for j=1:length(dataAvg(:,k));try dataAvg(j.*(pRel(j,k)<.5),k)=NaN;end;end;end
                     errorbar(ax7,dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2)),L,U,'color',['k'],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style);
                    set(ax7,'ylim',[-2,20]);
                end
                % Plot Phases
                hold on;
                if max(max(dataPhase))>2*pi
                    disp(['dataPhase > 2Pi, shifting by: ' num2str(-(max(max(dataPhase))-2*pi))]);
                    dataPhase=dataPhase-(max(max(dataPhase))-2*pi);
                end
                % can add 2Pi and maintain relationship with A injector
                if mean(mean(dataPhase)) <=-pi
                    dataPhase= dataPhase+2*pi;
                    disp('Shifting up by 2pi');
                end
                if mean(mean(dataPhase)) >=pi
                    dataPhase= dataPhase-2*pi;
                    disp('Shifting down by 2pi');
                end
                for i=1:size(doubleplot,1)
                    if Analysis==1
                        plot(dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,'-*','color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
                        ylabel('Phase [deg]');
                    elseif Analysis==2
                        for j=1:length(dataPhase(:,i));try dataPhase(j.*(pRel(j,i)<CutPow),i)=NaN;end;end
                        phaseH(:,i)=plot(ax8,dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,'-*','color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
                        ylabel(ax8,'[deg]');set(ax8,'ylim',[-400,400]);set(ax8,'xticklabel',[]);
                        set(ax8,'xlim',xlim);
                    end
                end
                %  Plot Fft Spectrum
                plot(ax9,dat(1).impacts(1:size(data,2)),100*pRel(:,1),'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
                plot(ax9,dat(1).impacts(1:size(data,2)),100*pRel(:,2),'-*','color', in(n).color{2}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
                 ylabel(ax9,'[%]');set(ax9,'ylim',[0,100]);
                 set(ax9,'xlim',xlim); title(ax9,'Inj. Mode % of Reconstruction');
                plot(ax9,xlim,[CutPow,CutPow].*100,'--k');
                figure(h2);
            elseif in(n).fftPlot & Analysis ==2
                % flow
                errorbar(ax7,dat(1).impacts(1:size(data,2)),param(:,2),param(:,3),param(:,3),'color',['k'],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style);
                set(ax7,'ylim',[-2,15]);
                
                % Plot Phases
                hold on;
                if max(max(dataPhase))>2*pi
                    dataPhase=dataPhase-(max(max(dataPhase))-2*pi);
                end
                for j=1:length(dataPhase(:,1));try dataPhase(j.*(pRel(j,i)<CutPow),1)=NaN;end;end
                plot(ax8,dat(1).impacts(1:size(data,2)),dataPhase(:,1).*180./pi,'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
                ylabel(ax8,'[deg]');set(ax8,'ylim',[0,400]);set(ax8,'xticklabel',[]);
                set(ax8,'xlim',xlim);
               
                 %  Plot Fft Spectrum
                plot(ax9,dat(1).impacts(1:size(data,2)),100*pRel(:,1),'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
                plot(ax9,dat(1).impacts(1:size(data,2)),100*pRel(:,2),'-*','color', in(n).color{2}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
                 ylabel(ax9,'[%]');set(ax9,'ylim',[0,100]);
                 set(ax9,'xlim',xlim); title(ax9,'Inj. Mode % of Reconstruction');
                plot(ax9,xlim,[CutPow,CutPow].*100,'--k');
                figure(h2);
            end
            %plot(dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2)),'k','LineWidth', lnwdth, 'LineStyle', in(n).style);
            if Analysis==1 & in(n).doubleplot==1
                plot(dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2)),'color' ,'k','marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style);
                plot(xlim,[0,0],'--k')
                ylabel('Toroidal Flow [km/s]'); set(gca,'ylim',[-10,10]);
                xlabel('Impacts [cm]');
            elseif Analysis==2 %& in(n).doubleplot==1
               % plot(ax7,dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2)),'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style);
                plot(ax7,xlim,[0,0],'--k'); set(ax7,'xlim',xlim);
                ylabel(ax7,'[km/s]'); set(ax7,'ylim',[-11,11]);set(ax7,'xticklabel',[]);
                h2.delete;
                linkaxes([ax6 ax7 ax8 ],'x');
            end

            figure(h) % make other current
        elseif plotType ==2
            % find cycle averaged temperature and fluxuation amplitude
            
        
        
        end
    end
    
    %% offset each line for plot 1
    for j = 1:size(data, 2)
       %PhaseVelocity.Velocity(j,:) = PhaseVelocity.Velocity(j,:)+ (j-1)*50;   
        data(:, j) = data(:, j) + (j-1) * offset;
        zeroline(:,j) = zeros(size(data,1)/2,1)+(j-1) * offset;
    end

    %% Plot Data
    if in(n).error
        time = ndgrid(dat(1).time, 1:size(data, 2));
        t(n, :) = errorbar(ax, time, data, errorL, errorU, 'color', in(n).color, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
%         for m = 1:size(t, 2)
%             errorbar_tick(t(n, m), errWdth); % adjust errorbar width
%         end
    else
        if ~isempty(in(n).doubleplot) && ~plotTor % plot both fibers
            t(n, :) = plot(ax, time(1:length(dat(1).time)), data(1:length(dat(1).time),:), ...
                'color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
            t(n, :) = plot(ax, time(length(dat(1).time)+1:end)', data(length(dat(1).time)+1:end,:), ...
                'color', in(n).color{2}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
        elseif  ~isempty(in(n).doubleplot) && plotTor % plot the toroidal flow
            t(n, :) = plot(ax, time(1:length(dat(1).time)), data(length(dat(1).time)+1:end,:)-data(1:length(dat(1).time),:) +zeroline, ...
                'color', 'k', 'LineWidth', lnwdth, 'LineStyle', in(n).style);
            plot(ax, [time(1),time(length(dat(1).time))], zeroline([1,size(data,1)./2],:), ...
                '--k', 'LineWidth', .5, 'LineStyle', in(n).style);
        else
            t(n, :) = plot(ax, time, data, 'color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
        end
    end
    
end

figure(h); % make first figure current

%% Impact Parameter Labels on Right
for n = 1:size(data, 2)
    y = offset * (n-1) + 0.1 * offset;
    text(timebound(end) + 0.02 * (timebound(end) - timebound(1)), y, num2str(dat(1).impacts(n), 2), 'fontsize', fntsz);
    plot(ax, [dat(1).time(1), dat(1).time(end)], zeros(2) + (n-1) * offset, '-', 'color', 'k');
end

% %% Legend
% legendText = cell(1, length(in)); % initialize
% legendHands = zeros(1, length(in)); % initialize
% for n = 1:length(in)
%     legendText{n} = in(n).legend;
%     legendHands(n) = t(n, 1);
% end
% legend(ax, legendHands, 'Location', 'NorthEastOutside',...
%      legendText, 'fontsize', fntsz);

 %% Misc. Figure Properties
if ~compactCurrents
    xlabel('Time [ms]','fontsize', fntsz);
else
    set(gca,'xticklabel',[]);
end
ylabel(sidebar, 'fontsize', fntsz);

title([in(1).legend  ':' titles], 'fontsize', fntsz);
%set(ax, 'XLim', timebound);
set(gca, 'LineWidth', lnwdth);
set(gca, 'fontsize', fntsz);
box on;
grid on;
if or(plotType == 2, plotType == 3) % temperature or Intensity
    yLowerLim = 0;
    yUpperLim = offset * size(data, 2);
else
    yLowerLim = -0.5 * offset;
    yUpperLim = offset * size(data, 2) - 0.5 * offset;
end
set(gca, 'YLim', [yLowerLim, yUpperLim]);
set(gca, 'YTick', []);

nt = text(timebound(end) + 0.125 * (timebound(end) - timebound(1)), (offset * size(data, 2))/2, 'R  [cm]', 'fontsize', fntsz);
set(nt, 'rotation', -90)

%% Figure Properties for Averages
if plotAverages
    figure(h2) % make current
    
    title(['Mean ' titles ' and Fluctuations'], 'fontsize', fntsz);
    set(ax2, 'XLim', [dat(1).impacts(1) - 1, dat(1).impacts(end) + 1]);
    set(gca, 'LineWidth', lnwdth);
    set(gca, 'fontsize', fntsz);
    box on;
    grid on;
    xlabel('Major Radius [cm]');
    ylabel(units);
    
    % Legend
    legendText = cell(1, length(in)); % initialize
    legendHands = zeros(1, length(in)); % initialize
    for n = 1:length(in)
        legendText{n} = in(n).legend;
        legendHands(n) = t2(n);
    end
    legend(ax2, legendHands, 'Location', 'NorthEastOutside',...
        legendText, 'fontsize', fntsz);
end

%% Plot Currents
if plotCurrents
    % plot currents
    if ~compactCurrents
        j = figure('Visible', 'on', 'Name', 'MULTIPLOT-Currents', 'Position',...
            [5, 35, figureWidth, analysisHeight], 'Color', [1 1 1]);
        ax2 = axes('Parent', j, 'Position', [0.15, 0.08, 0.8, 0.15], 'FontSize', fntsz+1); % currents
    else
        ax2 = axes('Parent', h, 'Position', [.075, .08, .8, .15], 'FontSize', fntsz+1);
    end
    if length(in) > 1
        load(['dat' num2str(in(1).shot) '10.mat']);
    end
    if dat(1).shotRef >999999
        plot(ax2, dat(1).iinjaTime.*in(1).injTimeScale, ...
            dat(1).iinja.*in(1).injScale,'LineWidth', lnwdth); 
        hold on;
        plot(ax2, dat(1).iinjbTime.*in(1).injTimeScale, ...
            dat(1).iinjb.*in(1).injScale, 'color', 'red', 'LineWidth', lnwdth);
        hold on;
        plot(ax2, dat(1).iinjcTime.*in(1).injTimeScale, ...
            dat(1).iinjc.*in(1).injScale, 'color', 'green', 'LineWidth', lnwdth);
        hold on;
    else
         plot(ax2, dat(1).iinjxTime.*in(1).injTimeScale, ...
            dat(1).iinjx.*in(1).injScale,'LineWidth', lnwdth); 
        hold on;
        plot(ax2, dat(1).iinjyTime.*in(1).injTimeScale, ...
            dat(1).iinjx.*in(1).injScale, 'color', 'blue', 'LineWidth', lnwdth);
        hold on
    end
    for i = 1:length(in)
        plot(ax2, dat(1).ItorTime.*in(1).injTimeScale, ...
            Itor(:,i).*in(1).injScale, 'color', in(i).color{1}, 'LineWidth', lnwdth);
    end
    set(gca, 'LineWidth', lnwdth);
    xlabel('Time [ms]');
    ylabel('I_{INJ} [kA]');
    linkaxes([ax,ax2],'x');
    set(ax2, 'XLim', timebound);
    grid on;
    
end
if in(1).shot == 151217023 || in(1).shot == 151217024 || in(1).shot == 151217025 || in(1).shot == 151217026
    cd(['T:\IDS\Analysis Repository\' num2str(in(1).shot) ]);
else
    cd(['T:\IDS\Analysis Repository\Alternate Phasings\60D2\']);
end
%saving
if saving
    pause(1.5);% necessary to get figure size correct
    saveas(h, ['Lines L' num2str(in(1).line) titles num2str(in(1).shot)], 'fig');
    saveas(h6, ['Analysis L' num2str(in(1).line) num2str(in(1).shot)], 'fig');
    saveas(h, ['Lines L' num2str(in(1).line) titles num2str(in(1).shot)], 'bmp');
    saveas(h6, ['Analysis L' num2str(in(1).line) num2str(in(1).shot)], 'bmp');
    if in(1).fftPlot
        saveas(h7, ['Spectrum L' num2str(in(1).line) num2str(in(1).shot)], 'fig');
        saveas(h7, ['Spectrum L' num2str(in(1).line) num2str(in(1).shot)], 'bmp');
    end
end
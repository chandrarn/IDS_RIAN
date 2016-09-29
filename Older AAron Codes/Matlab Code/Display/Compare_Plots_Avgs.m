% Basic code for comparing plots of NIMROD output and IDS data
%
% Updated November 1st, 2013 to display only one plot with time
% horizontally and space vertically.

clear all; close all; clc;

addpath('S:\MC_IDS\Matlab Code\Data Repository'); % where corrected IDS data is kept
addpath('S:\General Matlab'); % where corrected IDS data is kept

%% Plot Settings

%% Plot_1 Settings: %%%%

shot = 12949910;
shotb = 12949610;

plotType = 3; % 3 = IDS velocity
              % 4 = IDS temperature
              % 5 = IDS residual
              % 6 = IDS amplitude
              % 7 = IDS intensity (area * cal factor)
              % 8 = IDS residual / area
              % 9 = IDS calculated "delta R plasma"
              % 10 = Upper temperature error bar
              % 11 = Lower temperature error bar
              % 12 = Upper velocity error bar
              % 13 = Lower velocity error bar
shiftTime = 0; % shift time axis for plot a [ms]
shiftVel = -7; % km/s velocty shift
torPlot = 1; % true for toroidal array (fibers 1 to 36)
% false for poloidal array (fibers 37 to 72)

%% Options

deltaR = 0; % calculate the plasma displacement assuming temperature is false
uniformTemp = 20; % [eV]
skinny = 1; % make much shorter plot for zoomed in time slice
spa045 = 0; % also plot amperican loop current from 45 degrees

saveFigure = 0; % save figure to file
fileName = 'S:\HIT_SI\APS\2013\Poster\test'; % file name for .png image

timeInMs = 1; % displays time in ms, otherwise time point number
useImpacts = 1; % plot the x scale in terms of impact parameter, else channel

screenArea = 0; % run velocity and temperature through a filter to discard points with low area
scArea = 100; % lower limit for Area

screenSNR = 0; % discard points with high residual / area ratio
scSNR = 0.45; % upper limit for residual / area

findTLim = 0; % print indices to the command line which correspond to times
tLimMS = [.1 2]; % time points in ms

%% General Settings: %%%%

tempLim = [0 100];
velLim = [-10 10];
resLim = [0 150];
ampLim = [0 900];
intLim = [0 500];
snrLim = [0 0.5];
deltaRLim = [0 60];
tempErrLim = [0 20];
velErrLim = [0 5];

jLim = [-35 100]; % Limit for current plots

timeLim = [1.4 2.0];

chan_ranget = [8:24]; % toroidal, mohawk port in midplane
% chan_ranget = [8:28]; % toroidal, mohawk port perp.
% chan_ranget = [7:29]; % toroidal, 71 degree port
% chan_ranget = [10:27]; % toroidal, axial port
% chan_ranget = 1:30; % NIMROD mohawk

chan_rangep = [46:63]; % poloidal

%% Calculate and Plot

eval(sprintf('load(''dat%i'');', shot)); % Real HIT-SI Data
shotRef = dat.shotRef;

if torPlot
    chan_range = chan_ranget;
    ylabi = dat.label1;
else
    chan_range = chan_rangep;
    ylabi = dat.label2;
end
% Trim all Data for channel range
dat = trimRange(dat, chan_range);
dat.time = dat.time + shiftTime;
%             dat.ItorTime = dat.ItorTime + shiftTime;

if timeInMs
    x = dat.time;
    xlab = {'Time [ms]'};
else
    x = 1:length(dat.time);
    xlab = {'Time Point Number'};
end
if useImpacts
    y = dat.impacts;
    ylab = ylabi;
else
    y = dat.peaks;
    ylab = 'Channel Number';
end
ylim = sort([y(1), y(end)]); % 'sort' ensures the limits are in ascending order

t_avg = 1e3 * mean(diff(dat.time));
title2 = ['Interval = ' num2str(t_avg, 3) ' \mu' 's'];
title3 = [', ' dat.title];
if deltaR
    [dat] = calcDeltaR(dat, uniformTemp, t_avg);
end
if screenArea
    [dat] = screenAreaFun(dat, scArea);
end
if screenSNR
    [dat] = screenSNRFun(dat, scSNR);
end
if findTLim % special feature for convertime time to index
    [~, n1] = min(abs(dat.time - tLimMS(1)));
    [~, n2] = min(abs(dat.time - tLimMS(2)));
    n1
    n2
end


switch plotType
    case 3
        % IDS Velocity
        data = shiftVel + dat.vel;
        title1 = 'Velocity [km/s]';
        clim = velLim;
        
    case 4
        % IDS Temperature
        data = dat.temp;
        title1 = 'Temperature [eV]';
        clim = tempLim;
        
    case 5
        % IDS Residual
        data = dat.residual;
        title1 = 'Residual [counts]';
        clim = resLim;
        
    case 6
        % IDS Amplitude
        data = squeeze(dat.fit_par(:, 1, :) ./ (sqrt(2*pi)*dat.fit_par(:, 3, :)));
        title1 = 'Amplitude [counts]';
        clim = ampLim;
        
    case 7
        % IDS Intensity
        data = dat.int;
        title1 = 'Intensity [Arb. U.]';
        clim = intLim;
        
    case 8
        % IDS Residual / Area
        data = dat.residual ./ squeeze(dat.fit_par(:, 1, :));
        title1 = 'Residual / Area [Arb. U.]';
        clim = snrLim;
        
end

%% Averaging

nTimeLim(1) = find(dat.time > timeLim(1), 1);
nTimeLim(2) = find(dat.time > timeLim(2), 1);

for n = 1:size(data, 2);
    dataAvg(n) = mean(data(~isnan(data(:, n)), n));
    dataStd(n) = std(data(~isnan(data(:, n)), n));
end

%% Shot B

clear dat;
eval(sprintf('load(''dat%i'');', shotb)); % Real HIT-SI Data b
shotRefB = dat.shotRef;

% Trim all Data for channel range
dat = trimRange(dat, chan_range);
dat.time = dat.time + shiftTime;

switch plotType
        
    case 3
        % IDS Velocity
        data = shiftVel + dat.vel;
        title1 = 'Velocity [km/s]';
        clim = velLim;
        
    case 4
        % IDS Temperature
        data = dat.temp;
        title1 = 'Temperature [eV]';
        clim = tempLim;
        
    case 5
        % IDS Residual
        data = dat.residual;
        title1 = 'Residual [counts]';
        clim = resLim;
        
    case 6
        % IDS Amplitude
        data = squeeze(dat.fit_par(:, 1, :) ./ (sqrt(2*pi)*dat.fit_par(:, 3, :)));
        title1 = 'Amplitude [counts]';
        clim = ampLim;
        
    case 7
        % IDS Intensity
        data = dat.int;
        title1 = 'Intensity [Arb. U.]';
        clim = intLim;
        
    case 8
        % IDS Residual / Area
        data = dat.residual ./ squeeze(dat.fit_par(:, 1, :));
        title1 = 'Residual / Area [Arb. U.]';
        clim = snrLim;
        
end

%% Averaging

nTimeLim(1) = find(dat.time > timeLim(1), 1);
nTimeLim(2) = find(dat.time > timeLim(2), 1);

for n = 1:size(data, 2);
    dataAvgB(n) = mean(data(~isnan(data(:, n)), n));
    dataStdB(n) = std(data(~isnan(data(:, n)), n));
end

%% Plotting
fntsz = 24; % Font Size
lnwdth = 2;

n_skip = 3; % number of points to skip over

S = get(0,'ScreenSize');

h1 = figure('Visible','on','Name','Data Comparison','Position',[.3*S(3) 35 .4*S(3) S(4)-110],...
    'Color', [1 1 1]);
ax(1) = axes('Parent', h1, 'Position', [.2 .3 .6 .63], 'FontSize', fntsz);

eps = 0.3; % offset of plots

p(1) = plot(dataAvg, y, 'd-r');
hold on;
p(2) = plot(dataAvgB, y+eps, 'd-b');
set(p, 'LineWidth', lnwdth);

for n = 1:size(data, 2)
    p(3) = plot([dataAvg(n) - dataStd(n), dataAvg(n) + dataStd(n)], [y(n), y(n)], '-r');
    p(4) = plot([dataAvgB(n) - dataStdB(n), dataAvgB(n) + dataStdB(n)], [y(n)+eps, y(n)+eps], '-b');
    set(p, 'LineWidth', lnwdth);
end

grid on;
box on;
if skinny
    ylabel(ylab);
end
legend(p(1:2), num2str(shotRef), num2str(shotRefB), 'Location', 'SouthWest');
xlabel(title1);
set(gca, 'YLim', [ylim(1)-.5, ylim(2)+.5], 'LineWidth', lnwdth);
set(gca, 'XLim', clim);
% if timeInMs
%     set(gca, 'XLim', timeLim);
% end

% title([title3]);

%% Save

if saveFigure
    fig_save = getframe(h1);
    [Xfig, mapfig] = frame2im(fig_save);
    imwrite(Xfig, [fileName '.png']);
end








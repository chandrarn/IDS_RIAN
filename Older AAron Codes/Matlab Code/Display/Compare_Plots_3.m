% Basic code for comparing plots of NIMROD output and IDS data
%
% Updated November 1st, 2013 to display only one plot with time
% horizontally and space vertically.

clear all; 
%close all; 
%clc;
cd('T:\IDS\Data Repository'); % where corrected IDS data is kept My data
%addpath('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data\TEMP');% NSTX
addpath('S:\General Matlab');
 %cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data')

%% Plot Settings

%% Plot_1 Settings: %%%%
format long g
shot = 14070801810

plotType =3; % 3 = IDS velocity
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
              % 14 = Displacement
              % 15 = Cumulative intigration
shiftTime = 0; % shift time axis for plot a [ms]
shiftVel = -20; % km/s velocty shift
                % -20 for 129810 - 129820
                % -7 for 129499, etc.
torPlot = 1; % true for toroidal array (fibers 1 to 36)
% false for poloidal array (fibers 37 to 72)

%% Options

deltaR = 0; % calculate the plasma displacement assuming temperature is false
uniformTemp = 20; % [eV]
skinny = 0; % make much shorter plot for zoomed in time slice
spa045 = 0; % also plot amperican loop current from 45 degrees

saveFigure = 0; % save figure to file
fileName = 'S:\General Exam\Images\profile_129499'; % file name for .png image

timeInMs =1; % displays time in ms, otherwise time point number
useImpacts = 1; % plot the x scale in terms of impact parameter, else channel

screenArea = 0; % run velocity and temperature through a filter to discard points with low area
scArea = 100; % lower limit for Area

screenSNR = 0; % discard points with high residual / area ratio
scSNR = 0.45; % upper limit for residual / area

findTLim = 0; % print indices to the command line which correspond to times
tLimMS = [.1 2]; % time points in ms

interval = 6.9; % c17ne frame interval, in uS

hitsi3 = 1; % hitsi3 data, with 3 injectors instead of 2

%% General Settings: %%%%

tempLim = [0 100];
velLim = [-10 20];
distLim = [-15 15];
resLim = [0 150];
ampLim = [0 900];
intLim = [0 200];
snrLim = [0 0.5];
deltaRLim = [0 60];
tempErrLim = [0 20];
velErrLim = [0 5];
cumLim = [-1 1]%[-40 40]%[-40 20]%[-500 1500];
cumTime = [173 223]; % trim time

jLim = [-50 100]; % Limit for current plots

timeLim = [0 2.5]%[1.0 1.9];
%timeLim = [1.3,2.1];
%timeLim= [];

chan_ranget = [8:24]; % toroidal, mohawk port in midplane
% chan_ranget = [8:28]; % toroidal, mohawk port perp.
% chan_ranget = [8:27]; % toroidal, 71 degree port
% chan_ranget = [10:27]; % toroidal, axial port
% chan_ranget = 1:30; % NIMROD mohawk

chan_rangep = [46:63]; % poloidal

%% Calculate and Plot
format long g
eval(sprintf('load(''dat%0.0f'');', shot)); % Real HIT-SI Data

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
    x = dat.time;%.*1e3;% TEMP MULTIPLY, SOME SHOTS IN S NOT mS
    xlab = {'Time [ms]'};
else 
    x = 1:length(dat.time);
    xlab = {'Time Point Number'};
end

if plotType == 15
    x = dat.time(cumTime(1):cumTime(2));
end

if useImpacts
    y = dat.impacts;
    ylab = ylabi;
else
    y = dat.peaks;
    ylab = 'Channel Number';
end
ylim = sort([y(1), y(end)]); % 'sort' ensures the limits are in ascending order

[X, Y] = meshgrid(x, y);

t_avg = 1e3 * mean(diff(dat.time),1);
title2 = ['\Delta' 't = ' num2str(t_avg, 3) ' \mu' 's'];
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
    case 1
        
    case 2
        
    case 3
        % IDS Velocity
        %TEMPORARY FIX
        %dat.vel=dat.vel(1:201,:);
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
        
    case 9
        % Calculated delta R plasma
        data = dat.deltaR;
        title1 = ['Calculated \delta' 'r_{plasma} [cm]'];
        clim = deltaRLim;
        
    case 10
        % IDS Temperature Upper Error Bar
        data = dat.tempU;
        title1 = 'Temperature Upper Error [eV]';
        clim = tempErrLim;
        
    case 11
        % IDS Temperature Lower Error Bar
        data = dat.tempL;
        title1 = 'Temperature Lower Error [eV]';
        clim = tempErrLim;
        
    case 12
        % IDS Velocity Upper Error Bar
        data = dat.velU;
        title1 = 'Velocity Upper Error [km/s]';
        clim = velErrLim;
        
    case 13
        % IDS Velocity Lower Error Bar
        data = dat.velL;
        title1 = 'Velocity Lower Error [km/s]';
        clim = velErrLim;
    case 14
        %total displacement: interval * velocity = displacement
        data = (dat.vel+shiftVel)*(interval*.1);% interval*1e-6  uS to S, /1e-5 km to cm
        title1 = ['Displacement per \Delta' 't [cm]'];
        clim = distLim;
    case 15
        data = dat.vel(cumTime(1):cumTime(2),:);
        %data(251,:)=mean([data(250,:);data(252,:)]);
        %data(142,:)=mean([data(141,:);data(143,:)]);
        %data( ~any(data,2), : ) = [];  %rows
        DatMean = mean(data)';
        %p = polyfit(1:length(data(:,1)),(sum(data(:,:),2)./size(dat.vel,2))',1);
        
        for(i=1:size(data,2))
            p = polyfit(1:length(data(:,i)),(data(:,i))',1);
            data(:,i)=data(:,i)-(p(1).*(1:length(data(:,i)))+p(2))';
            %data(:,i)=data(:,i)-DatMean(i);
        end
        %x= dat.time(cumTime(1):cumTime(2));
        data = cumtrapz(x,data).*.1;
        clim = cumLim;
        data= data.*100; % to convert m to cm
        xlab = {'Time [ms]'};
        title1 = ['Total displacement [cm]'];
end

%% Plotting
fntsz = 24; % Font Size
lnwdth = 2;

n_skip = 3; % number of points to skip over

S = get(0, 'ScreenSize');

if ~skinny
    h1 = figure('Visible','on','Name','Data Comparison','Position',[5 35 S(3)-12 S(4)-110],...
        'Color', [1 1 1]);
    ax(1) = axes('Parent', h1, 'Position', [.1 .3 .82 .63], 'FontSize', fntsz);
else
    h1 = figure('Visible','on','Name','Data Comparison','Position',[.25*S(3) 35 .5*S(3) S(4)-110],...
        'Color', [1 1 1]);
    ax(1) = axes('Parent', h1, 'Position', [.1 .3 .863 .63], 'FontSize', fntsz);
end

hold on;
size(X)
size(Y)
size(data')
surf(X, Y, data');
shading interp;
colormap jet;
if ~skinny
    colorbar('FontSize', fntsz);
end
grid on;
box on;
if ~skinny
    ylabel(ylab);
end
set(gca, 'YLim', ylim, 'XTickLabel', [], 'LineWidth', lnwdth);
if timeInMs
    set(gca, 'XLim', timeLim);
end

title([title1 ',' title2 title3]);
caxis(clim);

%% Amperian Current

if spa045
    HitTree=Tree('hitsi',dat.shotRef);
    itor045=NATIVEvalue(HitTree.getNode('\I_TOR_SPA045').getData().data.getDoubleArray);
    itor045Time= NATIVEvalue(HitTree.getNode('\I_TOR_SPA045').getDimensionAt(0).data());
    itor045Time = 1e3 * itor045Time; % convert to ms
    itor045 = 1e-3 * itor045; % convert to kA
end

%% Plot Currents

ax(2) = axes('Parent', h1, 'Position', [.1 .08 .863 .19], 'FontSize', fntsz);

plt1 = '-k';
plt2 = '-r';
plt3 = '-b';
plt4 = '-g';

ha = plot(dat.ItorTime + shiftTime, dat.Itor, plt1);
set(ha, 'LineWidth', lnwdth);
hold on;
box on;

%plot injectors, 2/3
if ~hitsi3 
    hb = plot(dat.iinjxTime + shiftTime, dat.iinjx, plt2);
    hc = plot(dat.iinjyTime + shiftTime, dat.iinjy, plt3);
    set(hb, 'LineWidth', lnwdth);
    set(hc, 'LineWidth', lnwdth);
else
    hb = plot(dat.iinjaTime + shiftTime, dat.iinja, plt2);
    hc = plot(dat.iinjbTime + shiftTime, dat.iinjb, plt3);
    hd = plot(dat.iinjcTime + shiftTime, dat.iinjc, plt4);
    set(hb, 'LineWidth', lnwdth);
    set(hc, 'LineWidth', lnwdth);
    set(hd, 'LineWidth', lnwdth);
end

if spa045
    he = plot(itor045Time, itor045, plt4);
    set(he, 'LineWidth', lnwdth);
end



set(gca, 'YLim', jLim, 'YTick', [-100 -50 0 50 100]);
% set(gca, 'YLim', [1.2*min([min(dat.Itor) min(dat.iinjx) min(dat.iinjy)]) ...
%     1.2*max([max(dat.Itor) max(dat.iinjx) max(dat.iinjy)])]);

grid on;
if ~skinny
    if ~spa045
        if ~hitsi3
            legend('I_{TOR}', 'I_{INJ, X}', 'I_{INJ, Y}', 'Location', 'EastOutside');
        else
            legend('I_{TOR}', 'I_{INJ, A}', 'I_{INJ, B}', 'I_{INJ, C}', 'Location', 'EastOutside');
        end
    else
        legend('I_{TOR}', 'I_{INJ, X}', 'I_{INJ, Y}', 'I_{\phi = 45^o}', 'Location', 'EastOutside');
    end
end
ylabel('[kA]');
xlabel(xlab);

set(gca, 'XLim', timeLim, 'LineWidth', lnwdth);
linkaxes(ax, 'x');

%% Save

if saveFigure
    fig_save = getframe(h1);
    [Xfig, mapfig] = frame2im(fig_save);
    imwrite(Xfig, [fileName '.png']);
end








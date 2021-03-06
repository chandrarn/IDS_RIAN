% Basic code for comparing plots of NIMROD output and IDS data
%
% Updated November 1st, 2013 to display only one plot with time
% horizontally and space vertically.

clear all; close all; clc;
try 
    
    addpath('/home/aaron/IDS/Matlab/');
end
try
    %AddAllThePaths;
    addpath('T:\IDS\Data Repository\');
end

%% Plot Settings

%% Plot_1 Settings: %%%%
format long g
%shot = 15102210810;
shot = 15121702610;
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
              % 14 = Displacement
              % 15 = Cumulative intigration
% shiftTime = 1.335; % PSI-TET, shift time axis for plot a [ms]
% shiftTime = 0.9377; % NIMROD, shift time axis for plot a [ms]
line = 2; % which line to plot, if multiple
shiftTime = 0; % IDS data, shift time axis for plot a [ms]
shiftVel = -10; % km/s velocty shift
                % -20 for 129810 - 129820
                % -7 for 129499, etc. !!! REVISED TO -14
shiftVel2 = 0;-4.0376; % special case for second fiber array
torPlot = 2; % true for toroidal array (fibers 1 to 36)
% false for poloidal array (fibers 37 to 72)

%% Options

deltaR = 0; % calculate the plasma displacement assuming temperature is false
uniformTemp = 20; % [eV]
skinny = 0; % make much shorter plot for zoomed in time slice
spa045 = 0; % also plot amperican loop current from 45 degrees
spa000 = 1; % plot spa045 and spa000 from PDC tree
square = 0;
double = 0; % plot both arrays simultaniously, overrides torPlot


saveFigure = 0; % save figure to file
fileName = '/home/aaron/Dropbox/Thesis/Latex/Images/'; % file name for .png image

timeInMs = 1; % displays time in ms, otherwise time point number
useImpacts = 1; % plot the x scale in terms of impact parameter, else channel

screenArea = 0; % run velocity and temperature through a filter to discard points with low area
scArea = 100; % lower limit for Area

screenSNR = 0; % discard points with high residual / area ratio
scSNR = 0.45; % upper limit for residual / area

findTLim = 1; % print indices to the command line which correspond to times
tLimMS = [.1 2]; % time points in ms

hitsi3 = 1; % hitsi3 data, with 3 injectors instead of 2

%% General Settings: %%%%

tempLim = [0 30];
velLim = [-5 5];
distLim = [-15 15];
resLim = [0 150];
ampLim = [0 900];
intLim = [0 5];
snrLim = [0 0.5];
deltaRLim = [0 60];
tempErrLim = [0 20];
velErrLim = [0 5];
cumLim = [-15 15];%[-40 40]%[-40 20]%[-500 1500];
cumTime = []; % trim time

jLim = [-30 10]; % Limit for current plots

% timeLim = [1.5 2.3];
timeLim = [.7 2.0];
% timeLim = [1.664, 1.872]; % aligned for IDS and both codes
% timeLim = [0.64 0.97]; % testing NIMROD

chan_ranget = [12:23]; % toroidal, mohawk port in midplane
% chan_ranget = [8:28]; % toroidal, mohawk port perp.
% chan_ranget = [8:27]; % toroidal, 71 degree port
% chan_ranget = [8:24]; % toroidal, axial port
% chan_ranget = 1:30; % NIMROD mohawk

% chan_rangep = [46:63]; % poloidal
chan_rangep = [48:60]; % poloidal

%% Calculate and Plot
format long g
try % speed up runtime if data is already loaded
    dat;
catch
    eval(sprintf('load(''dat%0.0f'');', shot)); % Real HIT-SI Data
end

if torPlot == 1
    chan_range = chan_ranget;
    ylabi = dat(1).label1;
elseif torPlot == 0
    chan_range = chan_rangep;
    ylabi = dat(1).label2;
else
    chan_range = [chan_ranget chan_rangep];
    ylabi = dat(1).label1;
end

% Trim all Data for channel range
cd('T:\RChandra\NewCodes\IDS_Display');
dat = trimRange(dat, chan_range,1,[]);
%dat.time = dat.time + shiftTime;
%             dat.ItorTime = dat.ItorTime + shiftTime;

if timeInMs
    dat(1).time = dat(1).time;
    dat(1).time = dat(1).time + shiftTime;
    x = dat(1).time(1:length(dat(1).vel)); % TEMP MULTIPLY, SOME SHOTS IN S NOT mS
    %x = 1:size(dat.vel,1);
    xlab = {'time [ms]'};
else 
    %x = 1:length(dat.time);
    x = 1:size(dat(1).vel,1);
    xlab = {'Time Point Number'};
end

if square && isempty(timeLim)
    timeLim = [dat(1).time(1) dat(1).time(end)];
end

if plotType == 15
    x = dat(1).time(cumTime(1):cumTime(2));
end

if double || torPlot == 2 % this takes care of the two fibers being flipped
    %dat(1).impacts(18)
    %dat(1).impacts(length(chan_ranget)+1:end) = dat(1).impacts(end:-1:length(chan_ranget)+1);
    %dat(1).impacts(18)
end
if useImpacts
    y = dat(1).impacts;
    ylab = ylabi;
else
    y = dat(1).peaks;
    ylab = 'Channel Number';
end
ylim = sort([y(1), y(end)])% 'sort' ensures the limits are in ascending order
%ylim = [.9 43.5];

[X, Y] = meshgrid(x, y);

t_avg =  mean(diff(dat(1).time)) * 1e3;
title2 = [' \Delta' 't = ' num2str(t_avg, 3) ' \mu' 's'];
if ~square
    title3 = [', ' dat.title];
else
    title3 = [];
end

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
    [~, n1] = min(abs(dat(1).time - tLimMS(1)));
    [~, n2] = min(abs(dat(1).time - tLimMS(2)));
    n1
    n2
end


switch plotType
    case 1
        
    case 2
        
    case 3
        % IDS Velocity
        %TEMPORARY FIX FOR 129817
       % dat.vel=dat.vel(1:201,:);
        data = shiftVel + dat(line).vel;
        title1 = 'Velocity [km/s]';
        clim = velLim;
        if torPlot
            data(:,length(chan_ranget)+1:end) = data(:,length(chan_ranget)+1:end) + shiftVel2;
        end
        
    case 4
        % IDS Temperature
        data = dat(line).temp;
        title1 = 'Temperature [eV]';
        clim = tempLim;
        
    case 5
        % IDS Residual
        data = dat(line).residual;
        title1 = 'Residual [counts]';
        clim = resLim;
        
    case 6
        % IDS Amplitude
        data = squeeze(dat(line).fit_par(:, 1, :) ./ (sqrt(2*pi)*dat(line).fit_par(:, 3, :)));
        title1 = 'Amplitude [counts]';
        clim = ampLim;
        
    case 7
        % IDS Intensity
        data = dat(line).int./100;
        title1 = 'Intensity [Arb. U.]';
        clim = intLim;
        
    case 8
        % IDS Residual / Area
        data = dat(line).residual ./ squeeze(dat(line).fit_par(:, 1, :));
        title1 = 'Residual / Area [Arb. U.]';
        clim = snrLim;
        
    case 9
        % Calculated delta R plasma
        data = dat(line).deltaR;
        title1 = ['Calculated \delta' 'r_{plasma} [cm]'];
        clim = deltaRLim;
        
    case 10
        % IDS Temperature Upper Error Bar
        data = dat(line).tempU;
        title1 = 'Temperature Upper Error [eV]';
        clim = tempErrLim;
        
    case 11
        % IDS Temperature Lower Error Bar
        data = dat(line).tempL;
        title1 = 'Temperature Lower Error [eV]';
        clim = tempErrLim;
        
    case 12
        % IDS Velocity Upper Error Bar
        data = real(dat(line).velU);
        title1 = 'Velocity Upper Error [km/s]';
        clim = velErrLim;
        
    case 13
        % IDS Velocity Lower Error Bar
        data = abs(dat(line).velL);
        title1 = 'Velocity Lower Error [km/s]';
        clim = velErrLim;
    case 14
        %total displacement: interval * velocity = displacement
        data = (dat(line).vel + shiftVel) * t_avg; % interval*1e-6  uS to S, /1e-5 km to cm
        title1 = ['Displacement per \Delta' 't [cm]'];
        clim = distLim;
    case 15
        data = dat(line).vel(cumTime(1):cumTime(2), :);
        %data(251,:)=mean([data(250,:);data(252,:)]);
        %data(142,:)=mean([data(141,:);data(143,:)]);
        %data( ~any(data,2), : ) = [];  %rows
        DatMean = mean(data)';
        p = polyfit(1:length(data(:, 1)), (sum(data(:, :), 2) ./ size(dat.vel, 2))', 1);
        
        for i = 1:size(data, 2)
            p = polyfit(1:length(data(:, i)), (data(:, i))', 1);
            data(:, i) = data(:, i) - (p(1) .* (1:length(data(:, i))) + p(2))';
            %data(:,i)=data(:,i)-DatMean(i);
        end
        x= dat(1).time(cumTime(1):cumTime(2));
        data = cumtrapz(x, data);
        clim = cumLim;
        data= data .* 100; % to convert m to cm
        xlab = {'Time [ms]'};
        title1 = 'Cumulative displacement [cm]';
end

%% Plotting
fntsz = 24; % Font Size
lnwdth = 1;

n_skip = 3; % number of points to skip over

S = get(0, 'ScreenSize');

if (~skinny & ~square & ~double)
    h1 = figure('Visible','on','Name','Data Comparison','Position',[5 35 S(3)-1000 S(4)-110],...
        'Color', [1 1 1]);
    ax(1) = axes('Parent', h1, 'Position', [.1 .3 .82 .63], 'FontSize', fntsz);
elseif square
    h1 = figure('Visible','on','Name','Data Comparison','Position', [.15*S(3) 0.1*S(4) .25*S(3) .7*S(4)],...
        'Color', [1 1 1]);
    ax(1) = axes('Parent', h1, 'Position', [0.08, 0.1, 0.98, 0.8], 'FontSize', fntsz);
elseif double
    h1 = figure('Visible','on','Name','Data Comparison','Position',[.25*S(3) 35 .5*S(3) S(4)-110],...
        'Color', [1 1 1]);
    ax(3) = axes('Parent', h1, 'Position', [.1 .3 .834 .315], 'FontSize', fntsz);
    ax(1) = axes('Parent', h1, 'Position', [.1 .65 .834 .315], 'FontSize', fntsz);
else
    h1 = figure('Visible','on','Name','Data Comparison','Position',[.25*S(3) 35 .5*S(3) S(4)-110],...
        'Color', [1 1 1]);
    ax(1) = axes('Parent', h1, 'Position', [.1 .3 .863 .63], 'FontSize', fntsz);
end

hold on;
%data=data(1:201,:);
if torPlot == 2 % special instructions to trim middle region
    data(:,length(chan_ranget)+1) = zeros(size(data,1),1)*NaN;
end
if double
    h2 = surf(ax(1),X(1:length(chan_ranget)+1,:), Y(1:length(chan_ranget)+1,:), data(:,1:length(chan_ranget)+1)');
    h3 = surf(ax(3),X(length(chan_ranget)+1:end,:), Y(length(chan_ranget)+1:end,:), data(:,length(chan_ranget)+1:end)');
    set(h2, 'LineWidth', lnwdth);
    set(h3, 'LineWidth', lnwdth);
    view(ax(1),[0 90]);
    view(ax(3),[ 0 90]);
    shading(ax(1),'interp');
    shading(ax(3),'interp');
    colormap(ax(1),'jet');
    colormap(ax(3),'jet');
    grid(ax(1),'on');
    grid(ax(3),'on');
    box(ax(1),'on');
    box(ax(3),'on');
    cb = colorbar('peer',ax(1),'FontSize', fntsz);
    cb1 = colorbar('peer',ax(3),'FontSize', fntsz);
    ylabel(ylab);
    caxis(ax(1),clim);
    caxis(ax(3),clim);
    set(ax(1), 'YLim', [5,40], 'LineWidth', lnwdth);
    set(ax(3), 'YLim', [-40,-15], 'LineWidth', lnwdth);
    xlabel(ax(3),'time [ms]');
    if timeInMs
        set(ax(1), 'XLim', timeLim);
        set(ax(3), 'XLim', timeLim);
    else
        set(gca,'XLim',[1 length(dat(1).time)]);
    end
    %set(ax(1),'xtick',[]);
    set(gca,'xtick',[0:.1:2])
else
    
h2 = surf(ax(1),X,Y,data');
shading interp;
colormap jet;
grid on;
box on;

if ~skinny
    cb = colorbar('FontSize', fntsz);
end

if ~skinny
    ylabel(ylab);
end
if square
    set(cb, 'Position', [0.9, 0.1, 0.03, 0.8]);
end
set(gca, 'YLim', ylim, 'LineWidth', lnwdth);
if ~square
    set(gca,'XTickLabel', []);
else
    xlabel('time [ms]');
end

if timeInMs
    set(gca, 'XLim', timeLim);
else
    set(gca,'XLim',[1 length(dat(1).time)]);
end

%title([title1 ',' title2 title3]);
%title(['O II, 465 nm' ',' title2 title3]);
title(['O II, 465 nm' ',' title2 ', shot 151022108']);
caxis(clim);
set(cb,'Ytick',[-5:2.5:5])
set(gca,'xtick',[0:.2:2])
set(gca,'ytick',[-40:10:40])
box on;
breakyaxis([ -18,15]);
text(2.2,15,'Velocity [km/s]','Rotation',-90,'Fontsize',fntsz);

end

%% Amperian Current
import MDSplus.*
if spa045
    HitTree = Tree('hitsi',dat(1).shotRef);
    itor045 = NATIVEvalue(HitTree.getNode('\I_TOR_SPA045').getData().data.getDoubleArray);
    itor045Time = NATIVEvalue(HitTree.getNode('\I_TOR_SPA045').getDimensionAt(0).data());
    itor045Time = 1e3 * itor045Time; % convert to ms
    itor045 = 1e-3 * itor045; % convert to kA
elseif spa000
    Conn = Connection('landau.hit');
    Conn.openTree('pdc3',dat(1).shotRef);
    itor045 = NATIVEvalue(Conn.get('\I_TOR_SPA045')) .*1e-3;
    itor045Time = (NATIVEvalue(Conn.get('dim_of(\I_TOR_SPA045)'))).*1e+3;
    itor000 = NATIVEvalue(Conn.get('\I_TOR_SPA000')).*1e-3;
    itor000Time = (NATIVEvalue(Conn.get('dim_of(\I_TOR_SPA000)'))).*1e+3;
end

%% Plot Currents
if ~square
    ax(2) = axes('Parent', h1, 'Position', [.1 .08 (S(3)-1033)./1000 .19], 'FontSize', fntsz);

    plt1 = '-k';
    plt2 = '-r';
    plt3 = '-b';
    plt4 = '-g';

    %ha = plot(dat(1).ItorTime + shiftTime, dat(1).Itor.*1e-3, plt1);
    %set(ha, 'LineWidth', lnwdth);
    hold on;
    box on;

    %plot injectors, 2/3
    if ~hitsi3 
        hb = plot(dat(1).iinjxTime + shiftTime, dat(1).iinjx, plt2);
        hc = plot(dat(1).iinjyTime + shiftTime, dat(1).iinjy, plt3);
        set(hb, 'LineWidth', lnwdth);
        set(hc, 'LineWidth', lnwdth);
    elseif ~spa000 % TEMPORARILY MULTIPLYING BY 1e-3*, SOME SHOTS DIDNT DO THIS IN BATCH CORRECT
        hb = plot(dat(1).iinjaTime + shiftTime, dat(1).iinja, plt2);
        hc = plot(dat(1).iinjbTime + shiftTime, dat(1).iinjb, plt3);
        hd = plot(dat(1).iinjcTime + shiftTime, dat(1).iinjc, plt4);
        set(hb, 'LineWidth', lnwdth);
        set(hc, 'LineWidth', lnwdth);
        set(hd, 'LineWidth', lnwdth);
    else
        hb = plot(dat(1).iinjaTime + shiftTime, dat(1).iinja, plt2);
%         set(hb, 'LineWidth', lnwdth);
    end

    if spa045
        he = plot(itor045Time, itor045, plt4);
        set(he, 'LineWidth', lnwdth);
    elseif spa000
        he = plot(itor000Time, itor000, plt3);
%         set(he, 'LineWidth', lnwdth);
        hf = plot(itor045Time, itor045, plt4);
%         set(hf, 'LineWidth', lnwdth);
    end



    set(gca, 'YLim', jLim, 'YTick', [-20 -10 0 10 20]);
    % set(gca, 'YLim', [1.2*min([min(dat.Itor) min(dat.iinjx) min(dat.iinjy)]) ...
    %     1.2*max([max(dat.Itor) max(dat.iinjx) max(dat.iinjy)])]);

    grid on;
    if ~skinny
        if ~spa045 & ~spa000
            if ~hitsi3
                h_legend=legend('I_{TOR}', 'I_{INJ, X}', 'I_{INJ, Y}', 'Location', 'EastOutside');
            else
                h_legend=legend('I_{TOR}', 'I_{INJ, A}', 'I_{INJ, B}', 'I_{INJ, C}', 'Location', 'EastOutside');
            end
        elseif spa000
            h_legend=legend( 'I_{INJ, A}',  'I_{\phi = 0^o}', 'I_{\phi = 45^o}', 'Location', 'EastOutside');
        else
            h_legend=legend('I_{TOR}', 'I_{INJ, X}', 'I_{INJ, Y}', 'I_{\phi = 45^o}', 'Location', 'EastOutside');

        end
    end
    ylabel('[kA]');
    xlabel(xlab);
    set(h_legend,'FontSize',14);
    set(gca, 'XLim', timeLim, 'LineWidth', lnwdth,'FontSize',fntsz);
    linkaxes(ax, 'x');
    %set(gca,'xtick',[0:.1:2])
end
%% Save

if saveFigure
%     eval('print -dsvg /home/aaron/Dropbox/Thesis/Latex/Images/SIM_nim_int.svg');
    
%     fig_save = getframe(h1);
%     [Xfig, mapfig] = frame2im(fig_save);
%     imwrite(Xfig, [fileName 'Shot_' num2str(shot) '_type' num2str(plotType) '.png']);
end








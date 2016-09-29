% Basic code for comparing plots of NIMROD output and IDS data

clear all; close all; clc;

addpath('S:\MC_IDS\Matlab Code\Data Repository'); % where corrected IDS data is kept

%% Plot Settings

%% Plot_1 Settings: %%%%

shot_a = 12949910;

plotType_a = 4; % 3 = IDS velocity
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
shiftTime_a = 0; % shift time axis for plot a [ms]
shiftVel_a = 0; % km/s velocty shift
torPlot_a = 1; % true for toroidal array (fibers 1 to 36)
% false for poloidal array (fibers 37 to 72)

%% Plot_2 Settings: %%%%

shot_b = 12949910; 

plotType_b = 7; 

% shiftTime_b = 0.456; % shift time axis for plot b [ms]
shiftTime_b = 0;

shiftVel_b = 0; % km/s velocty shift
torPlot_b = 1; % true for toroidal array (fibers 1 to 36)
% false for poloidal array (fibers 37 to 72)

%% Options

deltaR = 0; % calculate the plasma displacement assuming temperature is false
uniformTemp = 20; % [eV]

mirror_a = 0; % flip plot a left to right
diffSettings = 0; % use different settings for each surface plot.

saveFigure = 0; % save figure to file
fileName = 'S:\MC_IDS\Analysis Plots\temp_int_129499'; % file name for .png image

timeInMs = 1; % displays time in ms, otherwise time point number
useImpacts = 1; % plot the x scale in terms of impact parameter, else channel

screenArea = 0; % run velocity and temperature through a filter to discard points with low area
scArea = 100; % lower limit for Area
scArea_b = 100;

screenSNR = 0; % discard points with high residual / area ratio
scSNR = 0.45; % upper limit for residual / area
scSNR_b = 0.45;

findTLim = 0; % print indices to the command line which correspond to times
tLimMS = [.1 2]; % time points in ms

%% General Settings: %%%%

tempLim = [0 100];
velLim = [-20 20];
resLim = [0 150];
ampLim = [0 900];
intLim = [0 500];
snrLim = [0 0.5];
deltaRLim = [0 60];
tempErrLim = [0 20];
velErrLim = [0 5];

% timeLim = [1.45 1.75];
timeLim = [1 2.2];

chan_ranget = [8:24]; % toroidal, mohawk port in midplane
% chan_ranget = [8:28]; % toroidal, mohawk port perp.
% chan_ranget = [7:29]; % toroidal, 71 degree port
% chan_ranget = [10:27]; % toroidal, axial port
% chan_ranget = 1:30; % NIMROD mohawk

chan_rangep = [46:63]; % poloidal

%% Calculations


%% Calculate and Plot
for n = 1:2
    if n == 1 % Plot 1
        plotType = plotType_a;
        shiftVel = shiftVel_a;
        shiftTime = shiftTime_a;
        shot = shot_a;
        torPlot = torPlot_a;
        

    else % Plot 2
        if diffSettings
            scSNR = scSNR_b;
            scArea = scArea_b;
        end

        plotType = plotType_b;
        shiftVel = shiftVel_b;
        shiftTime = shiftTime_b;
        shot = shot_b;
        torPlot = torPlot_b;

    end
    
    eval(sprintf('load(''dat%i'');', shot)); % Real HIT-SI Data
    
    if torPlot
        chan_range = chan_ranget;
        xlabi = dat.label1;
    else
        chan_range = chan_rangep;
        xlabi = dat.label2;
    end

%     switch plotType
%         case {1, 2}
%             nim.time = nim.time + shiftTime;
%             if timeInMs
%                 [X, Y] = meshgrid(chan_range(1):chan_range(2), nim.time);
%                 ylab = {'Time [ms]'};
%             else
%                 [X, Y] = meshgrid(chan_range(1):chan_range(2), 1:length(nim.time));
%                 ylab = {'Time Point Number'};
%             end
%             if NIMfSum == 0
%                 t2str = 'Time Step = ';
%             else
%                 t2str = 'Time Avg. = ';
%             end
%             title2 = [t2str num2str(nim.t_avg, 3) ' +/- ' num2str(nim.stdev, 2) ' \mu' 's'];
%             title3 = [];
% 
%         case {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}
            % Trim all Data for channel range
            dat = trimRange(dat, chan_range);
            dat.time = dat.time + shiftTime;
%             dat.ItorTime = dat.ItorTime + shiftTime;

            if timeInMs
                y = dat.time;
                ylab = {'Time [ms]'};
            else
                y = 1:length(dat.time);
                ylab = {'Time Point Number'};
            end
            if useImpacts
                x = dat.impacts;
                xlab = xlabi;
            else
                x = dat.peaks;
                xlab = 'Channel Number';
            end
            xlim = sort([x(1), x(end)]); % 'sort' ensures the limits are in ascending order

            [X, Y] = meshgrid(x, y);

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
%     end


    switch plotType
        case 1

        case 2

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
    end

%% Plotting
    fntsz = 22; % Font Size

    n_skip = 3; % number of points to skip over

    if n == 1
        S = get(0,'ScreenSize');
        h1 = figure('Visible','on','Name','Data Comparison','Position',[5 35 S(3)-12 S(4)-110],...
            'Color', [1 1 1]);
        ax(n) = axes('Parent', h1, 'Position', [.07 .12 .32 .79], 'FontSize', fntsz);
    else
        ax(n) = axes('Parent', h1, 'Position', [.45 .12 .32 .79], 'FontSize', fntsz);
        set(gca, 'YTickLabel', []);
    end

    hold on;
    surf(X, Y, data);
    shading interp;
    colormap jet;
    colorbar('FontSize', fntsz);
    grid on;
    box on;
    if n == 1
        ylabel(ylab);
    end
    %     set(gca, 'XTick', xtick, 'XTickLabel', xticklab);
    set(gca, 'XLim', xlim);
    if timeInMs
        set(gca, 'YLim', timeLim);
    end
    xlabel(xlab);
    title({[title1 ',']; [title2 title3]});
    caxis(clim);
    
    clear dat;

end

if torPlot_a == torPlot_b
    linkaxes(ax, 'xy');
else
    linkaxes(ax, 'y');
end

%% Plot Currents

ax(3) = axes('Parent', h1, 'Position', [.82 .12 .12 .79], 'FontSize', fntsz);

ha = [];
hb = [];
for n = 1:2
    if n == 1 % Plot 1
        shiftTime = shiftTime_a;
        shot = shot_a;
        plt = '-b';
        plt2 = '-c';
    else % Plot 2
        shiftTime = shiftTime_b;
        shot = shot_b;
        plt = '-r';
        plt2 = '-m';
    end
    
    eval(sprintf('load(''dat%i'');', shot)); % Real HIT-SI Data 
    ha(end+1) = plot(dat.Itor, dat.ItorTime + shiftTime, plt);
    hold on;
    box on;
    L{n} = dat.title;
    hb(end+1) = plot(dat.iinjx, dat.iinjxTime + shiftTime, plt);
    hb(end+1) = plot(dat.iinjy, dat.iinjyTime + shiftTime, plt2);
    % min1 = min(dat.Itor);
    % max1 = max(dat.Itor);

end

set(ha, 'LineWidth', 2);
set(hb, 'LineWidth', 2);
% 
% try
%     if L{1} == L{2}
%         legend(L{2});
%     end
% catch
%     legend(ha, L{1}, L{2});
% end


grid on;
set(gca, 'YTickLabel', []);
title('I_{tor} & I_{inj, x}');
xlabel('Current [kA]');
% xlabel('Time [ms]');
% ylabel('I_{TOR} [kA]');
% xlimit(1) = min([min1 min2]);
% xlimit(2) = max([max1 max2]);
set(gca, 'YLim', timeLim);
% set(gca, 'XLim', 1.2*xlimit);
linkaxes(ax, 'y');

%% Save

if saveFigure
    fig_save = getframe(h1);
    [Xfig, mapfig] = frame2im(fig_save);
    imwrite(Xfig, [fileName '.png']);
end








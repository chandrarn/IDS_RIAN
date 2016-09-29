clear all; close all; clc;
addpath('S:\MC_IDS\Matlab Code\NIMROD');
addpath('S:\MC_IDS\Matlab Code\');
addpath('S:\MC_IDS');

shot = [126814, 120911]; % 126690, 126871
NIMfSum = 19; % 5
spSum = 0; % look at chord sum % 1
tSum = 0; % look at time pattern

saveFigure = 0; % save figure to file
fileName = 'S:\HIT_SI\APS\2012\Poster\IDS'; % file name for .png image

plotType = 3; % 1 = NIMROD velocity
              % 2 = NIMROD false temperature
              % 3 = IDS velocity
              % 4 = IDS temperature
              % 5 = IDS residual
              % 6 = IDS amplitude
              % 7 = IDS area
              % 8 = IDS residual / area

timeShift = 1.17; % ms for second plot
              
xlim = [1.5 1.8];
velLim = [-105 20];
velSep = 20;
areaLim = [-8e3 1e3];
areaSep = 1e3;

cylim = [-20 40];

%% Load Currents
mdsopen('landau.hit::hitsi', shot(1));
iinjx = 1e-3*mdsvalue('\I_INJ_X');
tiinjx = 1e3*mdsvalue('dim_of(\I_INJ_X)');
iinjy = 1e-3*mdsvalue('\I_INJ_Y');
tiinjy = 1e3*mdsvalue('dim_of(\I_INJ_Y)');
itor1 = 1e-3*mdsvalue('\I_TOR_SPAAVG');
titor1 = 1e3*mdsvalue('dim_of(\I_TOR_SPAAVG)');
mdsclose();

if length(shot) > 1
    mdsopen('landau.hit::hitsi', shot(2));
    itor2 = 1e-3*mdsvalue('\I_TOR_SPAAVG');
    titor2 = 1e3*mdsvalue('dim_of(\I_TOR_SPAAVG)');
    mdsclose();
end

%% Load IDS Data
if and(spSum, ~tSum)
    eval(sprintf('load(''dat%is'');', shot(1))); % HIT-SI IDS Data
elseif and(~spSum, tSum)
    eval(sprintf('load(''dat%it'');', shot(1))); % HIT-SI IDS Data
elseif and(spSum, tSum)
    eval(sprintf('load(''dat%ist'');', shot(1))); % HIT-SI IDS Data
else
    eval(sprintf('load(''dat%i'');', shot(1))); % HIT-SI IDS Data
end

ch = dat.peaks;

%% Make Strings

txt = zeros(length(ch), 2);
load impacts;
for n = 1:length(ch)
    txt(n, :) = [impacts(dat.chSum(n, 1)), impacts(dat.chSum(n, end))];
end

fntsz = 22;
S = get(0, 'ScreenSize');
h1 = figure('Visible','on','Name','Data Comparison','Position',[5 6 S(3)-10 S(4)-75],...
    'Color', [1 1 1]);
ax(1) = axes('Parent', h1, 'Position', [.06 .38 .8 .56], 'FontSize', fntsz);

for m = 1:length(shot)
    if m == 2 % load second set of data
        if and(spSum, ~tSum)
            load(['nim' num2str(shot(2)) '_' num2str(NIMfSum) 's']); % NIMROD Calculation
            nim.time = nim.time + timeShift;
%             eval(sprintf('load(''dat%is'');', shot(m))); % HIT-SI IDS Data
%         elseif and(~spSum, tSum)
%             eval(sprintf('load(''dat%it'');', shot(m))); % HIT-SI IDS Data
%         elseif and(spSum, tSum)
%             eval(sprintf('load(''dat%ist'');', shot(m))); % HIT-SI IDS Data
%         else
%             eval(sprintf('load(''dat%i'');', shot(m))); % HIT-SI IDS Data
        end
        col = '-r';
    else
        col = '-b+';
    end

    switch plotType
        case 1
            % NIMROD velocity
            data = nim.vel;
            title1 = 'NIMROD Velocity [km/s]';
            clim = velLim;

        case 2
            % NIMROD temperature
            data = nim.temp;
            title1 = 'NIMROD Temperature [eV]';
            clim = tempLim;

        case 3
            % IDS Velocity
            if m == 1
                data = dat.vel;
                t1 = 'IDS';
                t2 = 'Velocity';
                unit = '[km/s]';
                ysep = velSep;
                clim = velLim;
            else
                data = nim.vel;
            end

        case 4
            % IDS Temperature
            data = dat.temp;
            title1 = 'IDS Temperature [eV]';
            clim = tempLim;

        case 5
            % IDS Residual
            data = dat.residual;
            title1 = 'IDS Residual [counts]';
            clim = resLim;

        case 6
            % IDS Amplitude
            data = squeeze(dat.fit_par(:, 1, :) ./ (sqrt(2*pi)*dat.fit_par(:, 3, :)));
            title1 = 'IDS Amplitude [counts]';
            clim = ampLim;

        case 7
            % IDS Area
            data = squeeze(dat.fit_par(:, 1, :));
            t1 = 'IDS';
            t2 = 'Fit Area';
            unit = '[Arb. U.]';
            clim = areaLim;
            ysep = areaSep;
            baseline = mean(mean(data(1:20, :)));
            data = data - baseline;

        case 8
            % IDS Residual / Area
            data = dat.residual ./ squeeze(dat.fit_par(:, 1, :));
            title1 = 'IDS Residual / Area [Arb. U.]';
            clim = snrLim;

    end
    
    for n = 1:length(ch)
        yticks(n) = -(n-1)*ysep;
    end
    yticks = yticks(end:-1:1);

%% Data Plotting
    if m == 1
        h = plot(dat.time, data(:, n), col);
    else
        h = plot(nim.time, data(:, n), col);
    end
    hold on;
    grid on;
    set(h, 'LineWidth', 2);
    text(xlim(2)*1.022, 0.9*ysep, 'Impact', 'FontSize', fntsz);
    text(xlim(2)*1.015, 0.6*ysep, 'Parameter', 'FontSize', fntsz);
    text(xlim(2)*1.01, 0, [num2str(txt(1,1), 2) ' to ' num2str(txt(1,2), 2) ' cm'],...
        'FontSize', fntsz);
    for n = 2:length(ch)
        if m == 1
            h = plot(dat.time, data(:, n) - (n-1)*ysep, col);
        else
            h = plot(nim.time, data(:, n) - (n-1)*ysep, col);
        end
        set(h, 'LineWidth', 2);
        text(xlim(2)*1.01, -(n-1)*ysep, [num2str(txt(n,1), 2) ' to ' num2str(txt(n,2), 2) ' cm'],...
            'FontSize', fntsz);
    end
end

set(gca, 'XLim', xlim, 'XTickLab', [], 'YTick', yticks, 'YTickLab', num2str(zeros(length(yticks), 1)),...
    'YLim', ylim);
ylabel({['Ion ' t2 ' ' unit]; ['Grid Spacing = ' num2str(ysep) ' ' unit]});
if length(shot) > 1
    title(['Spatially Averaged ' t1 ' ' t2 ', Shot ' num2str(shot(1)) ' and NIMROD']);
else
    title(['Spatially Averaged ' t1 ' ' t2 ', Shot ' num2str(shot)]);
end

%% Current Plot
ax(2) = axes('Parent', h1, 'Position', [.06 .12 .924 .2], 'FontSize', fntsz);
h = plot(tiinjx, iinjx, '-m');
hold on;
grid on;
set(h, 'LineWidth', 2);
if length(shot) == 1
    h = plot(tiinjy, iinjy, '-g');
else
    load NIM_IDS;
    nim.time = NIM_IDS.time + timeShift;
    h = plot(nim.time, NIM_IDS.iinjx, '-g');
end
set(h, 'LineWidth', 2);
h = plot(titor1, itor1, '-b');
set(h, 'LineWidth', 2);
if length(shot) > 1
%     load(['nim' num2str(shot(2)) '_0']); % NIMROD Calculation WITH UNALTERED TIME BASE
    h = plot(nim.time, nim.Itor, '-r');
    set(h, 'LineWidth', 2);
end
set(gca, 'XLim', xlim, 'YLim', cylim);
xlabel('Time [ms]');
ylabel('Current [kA]');
if length(shot) > 1
    legend('I_{INJ} X','I_{INJ} X NIM',['I_{TOR} ' num2str(shot(1))],['I_{TOR} NIMROD'],...
        'Location','EastOutside');
else
    legend('I_{INJ} X','I_{INJ} Y',['I_{TOR} ' num2str(shot(1))], 'Location', 'EastOutside');
end

linkaxes(ax, 'x');

%% Save

if saveFigure
    fig_save = getframe(h1);
    [Xfig, mapfig] = frame2im(fig_save);
    imwrite(Xfig, [fileName '.png']);
end

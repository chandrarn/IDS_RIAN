% Aaron Hossack
% Plot raw data and reconstructed fits

close all; clear all; clc;
addpath('T:/IDS/Data Analysis/');

%addAllThePaths;
%cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data')%NSTX data
cd('T:\IDS\Data Repository\');%my data
%% Settings

shot = 812951810;
chan_sets = [1:36; 37:72]; % available sets of channel numbers
chans1 = chan_sets(1, :);
deadCh = 5; % dead channel for displaying time instead of data

% Make Movie --------------------------------------------------------------
timeLim = [1:30];
saveMOVIE = 0;
%fileMOVIE = ['/home/aaron/IDS/Display/Images/fits' num2str(shot)];
fileMOVIE = strcat('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data\Temp\',num2str(shot));

timePt = 28; % time index
saveONEFRAME = 1;
fileONEFRAME = ['T:\IDS\Analysis Codes\Phase Data\Temp\' num2str(shot) 'time' num2str(timePt)];

%% Load and Prep Data

load(['dat' num2str(shot)]);

[n_time, n_pix, n_spatial] = size(dat.raw);

st = find(dat.peaks == chans1(1));
ed = find(dat.peaks == chans1(end));
chans = dat.peaks(st:ed); % exclude bad channels from channel range

%% MOVIE

fntsz = 10;
S = get(0,'ScreenSize');
h1 = figure('Visible','on','Name','IDS Fitting','Position',...
    [.005*S(3), .03*S(4), .99*S(3) .9*S(4)], 'Color', [1 1 1]);

n_chy = 4;
n_chx = 9;
g = reshape(chans1, n_chx, n_chy)'; % make grid of channel numbers corresponding to position in figure

k = 1;
for n = timeLim % Time Loop
    clf;
    
    for m = chans' % Channel Loop
        
        m1 = find(chans == m); % INDEX of channel number
        
        % Figure out where to put axes ------------------------------------
        
        [ny, nx] = find(g == m);
        ny = n_chy + 1 - ny; % invert 'ny' to count from the bottom
        
        sx = 1 / n_chx * (nx - 1) + .01; % lower left corner of x axis
        sy = 1 / n_chy * (ny - 1) + .007; % lower left corner of y axis
        
        dx = .88 * 1 / n_chx; % x width
        dy = .88 * 1 / n_chy; % y width
        
        % Surface Plots ---------------------------------------------------
        
        sy_1 = sy + 0.5 * dy; % start point is halfway up the allotted space for this channel
        dx_1 = 0.4 * dx; % half the width of the subplot
        dy_1 = 0.48 * dy; % half the width of the subplot
        
        ax(1) = axes('Parent', h1, 'Position', [sx sy_1 dx_1 dy_1], 'FontSize', fntsz);
        
        % Make Grid and plot Raw Data -------------------------------------
        
        xBound = dat.bounds(n, m1, 1) : dat.bounds(n, m1, 2);
        yBound = dat.bounds(n, m1, 3) : dat.bounds(n, m1, 4);
        
        [X, Y] = meshgrid(xBound, yBound);
        Z = squeeze(dat.raw(n, yBound, xBound));
        
        surf(X, Y, Z);
        shading interp;
        view([0 90]);
        hold on;
        offset = size(X, 2);
        set(gca, 'XLim', [xBound(1), xBound(1) + 3*offset], 'YLim', [yBound(1) yBound(end)]);
        
        % Make Grid for Fits and Guesses ----------------------------------
        
        xBoundf = linspace(dat.bounds(n, m1, 1), dat.bounds(n, m1, 2), 50);
        yBoundf = linspace(dat.bounds(n, m1, 3), dat.bounds(n, m1, 4), 150);

        [Xf, Yf] = meshgrid(xBoundf, yBoundf);
        
        % Reshape fine mesh for execution by Gaussian function
        
        xf(:, 1) = Xf(:);
        xf(:, 2) = Yf(:);
        
        % Evaluate Fits and Guesses ---------------------------------------
        
        zf = singletGauss2D(squeeze(dat.fit_par(n, m1, :)), xf); % calculate fit on fine mesh
        zg = singletGauss2D(squeeze(dat.guesses(n, m1, :)), xf); % calculate guess on fine mesh
        
        Zf = reshape(zf, size(Xf, 1), size(Xf, 2)); % reshape into 2D image
        Zg = reshape(zg, size(Xf, 1), size(Xf, 2)); % reshape into 2D image
        
        % Plot Fits and Guesses -------------------------------------------
        
        surf(offset + Xf, Yf, Zf); % Fit
        shading interp;
        surf(2 * offset + Xf, Yf, Zg); % Guess
        shading interp;
        set(gca, 'XTickLabel', [], 'YTickLabel', []);
        title(['Channel ' num2str(m)]);
        
        % Display Parameters ----------------------------------------------
        
        heights = linspace(.95, .05, 6);
        labels = {'Vol', 'x_0', 'y_0', '\sigma_x', '\sigma_y', 'Off'};
        
        for p = 1:6
            text('Units', 'normalized', 'Position', [1.05, heights(p)], 'String', ...
                labels{p}); % Labels
            text('Units', 'normalized', 'Position', [1.6, heights(p)], 'String', ...
                num2str(dat.guesses(n,m1,p),3)); % Initial Guesses
            text('Units', 'normalized', 'Position', [2.1, heights(p)], 'String', ...
                num2str(dat.fit_par(n,m1,p),3)); % Final Fits
        end
    
        % y Fit -----------------------------------------------------------
        
        ax(2) = axes('Parent', h1, 'Position', [sx sy dx dy_1], 'FontSize', fntsz);
        
        xind = ceil(length(xBound) / 2); % index
        xindf = ceil(length(xBoundf) / 2); % index
        
        plot(yBound, Z(:, xind), '+r');
        hold on;
        plot(yBoundf, Zf(:, xindf), '-b');
        plot(yBoundf, Zg(:, xindf), '-c');
        
        set(gca, 'XTickLabel', []);
        
        % Display Final Answers -------------------------------------------
        
        text('Units', 'normalized', 'Position', [0.02, 0.9], 'String', ...
            [num2str(dat.temp(n,m1),2) ' eV']); % Temperature
            
        text('Units', 'normalized', 'Position', [0.7, 0.9], 'String', ...
            [num2str(dat.vel(n,m1),2) ' km/s']); % Velocity
        
        if m == deadCh - 1

            % Display time in dead channel space --------------------------

            text('Parent', ax(1), 'Units', 'normalized', 'Position', [3.3, .6], 'String', ...
                ['Time Pt. ' num2str(n)]);
            text('Parent', ax(1), 'Units', 'normalized', 'Position', [3.2, .4], 'String', ...
                ['Time = ' num2str(dat.time(n), 4) ' ms']);
        end

    end % Channel Loop
    
    

    pause(0.1);

    % Save Frame ---------------------------
    gcf;
    F(k) = getframe(h1);

    if and(saveONEFRAME, n == timePt)
        [Xfig, mapfig] = frame2im(F(k));
        imwrite(Xfig, [fileONEFRAME '.png']);
    end

    k = k+1;
    clear ax;
    
end % Time Loop

if saveMOVIE
    % Save Movie
    movie2avi(F, [fileMOVIE '.avi'], 'FPS', 5, 'compression', 'none', 'quality', 75);
end
    

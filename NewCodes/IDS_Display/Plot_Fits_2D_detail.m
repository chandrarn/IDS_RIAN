% Aaron Hossack
% Plot raw data and reconstructed fits

%notes: if the channel numbers being used is odd, may need to change
%chan_sets, and n_chx,n_chy, (n_chx*n_chy = length(chans1)) 

%Cyan is initial guess, blue is final fit

%close all; 
clear all; clc;
addpath('~/IDS/Matlab/');
addpath('T:\IDS\Data Repository\');
%AddAllThePaths;

%% Settings

% shot = 12858010;
shot = 16101802210;
chan_sets = [3:32;32:61]; % available sets of channel numbers
chans1 = chan_sets(2, :);
deadCh = [4]; % dead channel for displaying time instead of data
line = 2; % line INDEX, ie: dat(#).vel

% Make Movie --------------------------------------------------------------
% timeLim = 30:48;
timeLim = [50:181];
saveMOVIE = 1;
fileMOVIE = ['T:\IDS\Data Repository\Fits ' num2str(shot)];

timePt = 100; % time index
saveONEFRAME = 0;
fileONEFRAME = ['T:\IDS\Analysis Repository\' num2str(shot) 'time' num2str(timePt) 'Chanset' num2str(chans1(1))];

%% Load and Prep Data
  
load(['dat' num2str(shot)]);

[n_time, n_pix, n_spatial] = size(dat(1).raw);
if isempty(timeLim)
    timeLim = 1:n_time;
end

st = find(dat(1).peaks == chans1(1));
ed = find(dat(1).peaks == chans1(end));
chans = dat(1).peaks(st:ed); % exclude bad channels from channel range

%% MOVIE

fntsz = 10;
S = get(0,'ScreenSize');
h1 = figure('Visible','on','Name','IDS Fitting','Position',...
    [0.005*S(3), 0.03*S(4), 0.99*S(3) 0.8*S(4)], 'Color', [1 1 1]);

n_chy = 5;
n_chx = 6;
g = reshape(chans1, n_chx, n_chy)'; % make grid of channel numbers corresponding to position in figure

k = 1;
for n = timeLim % Time Loop
    n
    clf;
    
    for m = chans' % Channel Loop
        
        m1 = find(dat(1).peaks == m); % INDEX of channel number
        
        % Figure out where to put axes ------------------------------------
        
        [ny, nx] = find(g == m);
        ny = n_chy + 1 - ny; % invert 'ny' to count from the bottom
        
        sx = 1 / n_chx * (nx - 1) + 0.01; % lower left corner of x axis
        sy = 1 / n_chy * (ny - 1) + 0.007; % lower left corner of y axis
        
        dx = 0.88 * 1 / n_chx; % x width
        dy = 0.88 * 1 / n_chy; % y width
        
        % Surface Plots ---------------------------------------------------
        
        sy_1 = sy + 0.5 * dy; % start point is halfway up the allotted space for this channel
        dx_1 = 0.4 * dx; % half the width of the subplot
        dy_1 = 0.48 * dy; % half the width of the subplot
        
        ax(1) = axes('Parent', h1, 'Position', [sx sy_1 dx_1 dy_1], 'FontSize', fntsz);
        
        % Make Grid and plot Raw Data -------------------------------------
        
        xBound = dat(line).bounds(n, m1, 1) : dat(line).bounds(n, m1, 2);
        yBound = dat(line).bounds(n, m1, 3) : dat(line).bounds(n, m1, 4);
        
        [X, Y] = meshgrid(xBound, yBound);
        Z = squeeze(dat(1).raw(n, yBound, xBound));
        
        surf(X, Y, Z);
        shading interp;
        view([0 90]);
        hold on;
        offset = size(X, 2);
        set(gca, 'XLim', [xBound(1), xBound(1) + 3*offset], 'YLim', [yBound(1) yBound(end)]);
        
        % Make Grid for Fits and Guesses ----------------------------------
        
        xBoundf = linspace(dat(line).bounds(n, m1, 1), dat(line).bounds(n, m1, 2), 50);
        yBoundf = linspace(dat(line).bounds(n, m1, 3), dat(line).bounds(n, m1, 4), 150);

        [Xf, Yf] = meshgrid(xBoundf, yBoundf);
        
        % Reshape fine mesh for execution by Gaussian function
        
        xf(:, 1) = Xf(:);
        xf(:, 2) = Yf(:);
        
        % Evaluate Fits and Guesses ---------------------------------------
        
        zf = singletGauss2D(squeeze(dat(line).fit_par(n, m1, :)), xf); % calculate fit on fine mesh
        zg = singletGauss2D(squeeze(dat(line).guesses(n, m1, :)), xf); % calculate guess on fine mesh
        
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
                num2str(dat(line).guesses(n,m1,p),3)); % Initial Guesses
            text('Units', 'normalized', 'Position', [2.1, heights(p)], 'String', ...
                num2str(dat(line).fit_par(n,m1,p),3)); % Final Fits
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
            [num2str(dat(line).temp(n,m1),2) ' eV']); % Temperature
            
        text('Units', 'normalized', 'Position', [0.7, 0.9], 'String', ...
            [num2str(dat(line).vel(n,m1),2) ' km/s']); % Velocity
        %display([num2str(n) ',' num2str(m1) ',' num2str(dat(line).vel(n,m1))]);
        
        if m == deadCh - 1

            % Display time in dead channel space --------------------------

            text('Parent', ax(1), 'Units', 'normalized', 'Position', [3.3, .6], 'String', ...
                ['Time Pt. ' num2str(n)]);
            text('Parent', ax(1), 'Units', 'normalized', 'Position', [3.2, .4], 'String', ...
                ['Time = ' num2str(dat(1).time(n), 4) ' ms']);
        end

    end % Channel Loop
    
    

    %pause(0.1);

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
    v=VideoWriter([fileMOVIE '.mj2'],'Motion JPEG 2000');
    %v.CompressionRatio=2;
    v.LosslessCompression=true;
%    v.Quality=100;
 %   v.VideoFormat='RGB24';
 v.FrameRate=2;
    
    open(v);
    writeVideo(v,F);
    close(v);
    %movie2avi(F, [fileMOVIE '.avi'], 'fps', 2, 'compression', 'none', 'quality', 100);
end
    

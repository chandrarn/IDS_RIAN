% Aaron Hossack
% First attempt at BD analysis of raw IDS data
%
% The '3' version does a better job- now, tight bounds are set around the
% toroidal and poloidal data, they are individually BD'd, and then placed
% onto an array of zeros.  Thus, the movie is the same size as the original
% so the channel numbers match.
%
% In this version (2), the BD analysis is done on the raw data, rather than
% the binned data.

clear all; close all; clc;

addpath('C:\Program Files\MDSplus\MATLAB');
addpath('S:\MC_IDS\Matlab Code\Core Codes\lsqcurvefit');
cd('S:\MC_IDS\Matlab Code\Core Codes');

%% Input Options -----------------
torPlot = 1; % 1 = save toroidal plots (if selected), 0 = save poloidal plots

% Weights -----------
plotWeights = 0;
saveWeights = 0;
fileWeights = 'S:\HIT_SI\Presentations\Friday Meetings\130503_IDS_BD_1\Weights_126690';

plotFraction = 0;
saveFraction = 0;
fileFraction = 'S:\HIT_SI\Presentations\Friday Meetings\130503_IDS_BD_1\Fraction';

% Topos -------------
plotTopo1 = 0;
saveTopo1 = 0;
fileTopo1 = 'S:\HIT_SI\Presentations\Friday Meetings\130503_IDS_BD_1\Topo1';

plotTopo2 = 0;
saveTopo2 = 0;
fileTopo2 = 'S:\HIT_SI\Presentations\Friday Meetings\130503_IDS_BD_1\Topo2';

plotTopo3 = 0;
saveTopo3 = 0;
fileTopo3 = 'S:\HIT_SI\Presentations\Friday Meetings\130503_IDS_BD_1\Topo3';

plotTopo4 = 0;
saveTopo4 = 0;
fileTopo4 = 'S:\HIT_SI\Presentations\Friday Meetings\130503_IDS_BD_1\Topo4';

plotTopo5 = 0;
saveTopo5 = 0;
fileTopo5 = 'S:\HIT_SI\Presentations\Friday Meetings\130503_IDS_BD_1\Topo5';

plotTopoComb = 0;
saveTopoComb = 0;
fileTopoComb = 'S:\HIT_SI\Presentations\Friday Meetings\130503_IDS_BD_1\TopoComb_126690';
combPixLim = [1 20]; % y limit for combination plot

% Chronos -----------
plotChrons = 0;
saveChrons = 0;
fileChrons = 'S:\HIT_SI\Presentations\Friday Meetings\130503_IDS_BD_1\Chrons_126690';

plotChronsW = 0;
saveChronsW = 0;
fileChronsW = 'S:\HIT_SI\Presentations\Friday Meetings\130503_IDS_BD_1\ChronsW_126690';

% Recombining Signals ------------
recombine = 0;
modes = 1:4; % array of mode numbers to recombine and add together
% falseShot = 13071902; % fake shot number, 'yymmdd##' format
% ^ Now outdated, switched false shot format

% General Settings ----------
shot = 129499;

% --- real space direction pixel bound - toroidal fiber
chan_bound_t = 30:180; % 129499
% chan_bound_t = 30:180; % 129591
% chan_bound_t = 60:170; % 129793, 129819

% --- real space direction pixel bound - poloidal fiber
chan_bound_p = 230:320; % 129499
% chan_bound_p = 210:320; % 129591
% chan_bound_p = 223:333; % 129793, 129819

% --- wavelength direction pixel bound
% lam_bound = 30:50; 
lam_bound = 29:49;

% --- array of time INDEX points for analysis
% time_bound = 70:320; % 129591
time_bound = 70:360; % 129499, 126690 w/ vac
% time_bound = 12:27; %
% time_bound = 10:29; %

falseShot = [num2str(shot) num2str(modes(end))]; % default - can be manually
% set to 'date' format if doing something other than a sum up to a certain
% mode number.

%%
% Load Constants and Parameters

t_override = 0;
% [param, options] = loadParams(shot, line);

% Load Data Array

cd('T:\PhantomMovies');
data = importdata(['shot' int2str(shot) '.mat']); % [counts] (time x wavelength space x channel space)
time = importdata(['t' int2str(shot) '.mat']); % [ms]
cd('S:\MC_IDS\Matlab Code\Core Codes');

% Unflip images to correct for flip during Python conversion

data = data(:, end:-1:1, end:-1:1);

% BD specific time bound
if isempty(time_bound)
    time_bound = 1:length(time);
end
data = data(time_bound, :, :);
time = time(time_bound);

% Record size for recombining
[n_time_o, n_pix_o, n_chan_o] = size(data);

% Account for the possibility of no wavelength bounds
if isempty(lam_bound)
    lam_bound = 1:n_pix_o;
end

% Account for empty channel bound
if isempty(chan_bound_t)
    chan_bound_t = 1:n_chan_o;
end
if isempty(chan_bound_p)
    chan_bound_p = 1:n_chan_o;
end

for m = 1:2
    if m == 1 % toroidal fiber
        chan_bound = chan_bound_t;
        if torPlot
            savePlot = 1;
        else
            savePlot = 0;
        end
    else % poloidal fiber
        chan_bound = chan_bound_p;
        if torPlot
            savePlot = 0;
        else
            savePlot = 1;
        end
    end
    
    % Crop Data Frame for BD
    data2 = data(:, lam_bound, chan_bound);
    
    % Arrange in vertical columns
    [n_time, n_pix, n_chan] = size(data2);
    BDdat = NaN*zeros(n_chan * n_pix, n_time);
    
    for n = 1:n_time
        BDdat(:, n) = reshape(squeeze(data2(n, :, :)), n_chan * n_pix, 1);
    end
    
    % Subtract off minimum value to reduce noise floor
    minVal = min(BDdat(:));
    
    disp(['Min = ' num2str(minVal)]);
    disp(['Median = ' num2str(median(BDdat(:)))]);
    disp(['Mean = ' num2str(mean(BDdat(:)))]);
    disp(['Std = ' num2str(std(BDdat(:)))]);
    
    BDdat = BDdat - minVal; % subtract minimum value in data
    
    %% Perform SVD -----------------------
    
    [U, S, V] = svd(BDdat, 'econ');
%     size(U)
%     size(S)
%     size(V)
    % Calculate Singular Values
    
    Ak = diag(S);
    
    % Calculate the cumulative amplitude squared fraction
    
    for n = 1:length(Ak)
        cumWeight(n) = sum(Ak(1:n).^2) / sum(Ak.^2);
    end
    
    % Rearrange data as it was before
    
    topos = zeros(n_time, n_pix, n_chan);
    
    for n = 1:n_time
        topos(n, :, :) = reshape(U(:, n), 1, n_pix, n_chan);
    end
    
    [X, Y] = meshgrid(1:n_chan, 1:n_pix);
    
    %% Plotting Settings
    
    S = get(0,'ScreenSize');
    lnwdth = 2;
    fntsz = 20;
    
    %% Plot Weights
    if plotWeights
        h1 = figure('Visible','on','Name','SVD Weights','Position',...
            [S(3)/4, S(4)/6, S(3)/2 2*S(4)/3], 'Color', [1 1 1]);
        h2 = axes('Parent', h1, 'FontSize', fntsz);
        h3 = semilogy(Ak, 'b-o');
        hold on;
        grid on;
        title(['SVD Weights, ''A_k'', shot ' num2str(shot)]);
        xlabel('Mode Number');
        ylabel('A_k');
        set(h3, 'LineWidth', lnwdth);
        set(gca, 'FontSize', fntsz);
        
        if saveWeights && savePlot
            fig_save = getframe(h1);
            [Xfig, mapfig] = frame2im(fig_save);
            imwrite(Xfig, [fileWeights '.png']);
        end
    end
    
    %% Plot Fraction
    if plotFraction
        h1 = figure('Visible','on','Name','Cumulative Fractional Weights','Position',...
            [S(3)/4, S(4)/6, S(3)/2 2*S(4)/3], 'Color', [1 1 1]);
        h2 = axes('Parent', h1, 'FontSize', fntsz);
        h3 = plot(cumWeight, 'b-o');
        hold on;
        grid on;
        title(['Cumulative Fractional Weight ' num2str(shot)]);
        xlabel('Mode Number');
        ylabel('A_{1 - k}^2 / A_{all}^2');
        set(h3, 'LineWidth', lnwdth);
        set(gca, 'FontSize', fntsz);
        
        if saveFraction && savePlot
            fig_save = getframe(h1);
            [Xfig, mapfig] = frame2im(fig_save);
            imwrite(Xfig, [fileFraction '.png']);
        end
    end
    
    %% Plot Topo 1
    if plotTopo1
        h1 = figure('Visible','on','Name','Topo 1','Position',...
            [S(3)/4, S(4)/6, S(3)/2 2*S(4)/3], 'Color', [1 1 1]);
        h2 = axes('Parent', h1, 'Position', [.1 .1 .8 .8], 'FontSize', fntsz);
        h3 = surf(X, Y, Ak(1)*squeeze(topos(1, :, :)));
        hold on;
        shading interp;
        colormap jet;
        colorbar;
        grid on;
        view([0 90]);
        title(['First Topo, shot ' num2str(shot)]);
        xlabel('Pixel Number (real space)');
        ylabel('Pixel Number (wavelength space)');
        set(h3, 'LineWidth', lnwdth);
        set(gca, 'FontSize', fntsz, 'XLim', [1 n_chan], 'YLim', [1 n_pix]);
        
        if saveTopo1 && savePlot
            fig_save = getframe(h1);
            [Xfig, mapfig] = frame2im(fig_save);
            imwrite(Xfig, [fileTopo1 '.png']);
        end
    end
    
    %% Plot Topo2
    if plotTopo2
        h1 = figure('Visible','on','Name','Topo 1','Position',...
            [S(3)/4, S(4)/6, S(3)/2 2*S(4)/3], 'Color', [1 1 1]);
        h2 = axes('Parent', h1, 'Position', [.1 .1 .8 .8], 'FontSize', fntsz);
        h3 = surf(X, Y, Ak(2)*squeeze(topos(2, :, :)));
        hold on;
        shading interp;
        colormap jet;
        colorbar;
        grid on;
        view([0 90]);
        title(['Second Topo, shot ' num2str(shot)]);
        xlabel('Pixel Number (real space)');
        ylabel('Pixel Number (wavelength space)');
        set(h3, 'LineWidth', lnwdth);
        set(gca, 'FontSize', fntsz, 'XLim', [1 n_chan], 'YLim', [1 n_pix]);
        
        if saveTopo2 && savePlot
            fig_save = getframe(h1);
            [Xfig, mapfig] = frame2im(fig_save);
            imwrite(Xfig, [fileTopo2 '.png']);
        end
    end
    
    %% Plot Topo3
    if plotTopo3
        h1 = figure('Visible','on','Name','Topo 3','Position',...
            [S(3)/4, S(4)/6, S(3)/2 2*S(4)/3], 'Color', [1 1 1]);
        h2 = axes('Parent', h1, 'Position', [.1 .1 .8 .8], 'FontSize', fntsz);
        h3 = surf(X, Y, Ak(3)*squeeze(topos(3, :, :)));
        hold on;
        shading interp;
        colormap jet;
        colorbar;
        grid on;
        view([0 90]);
        title(['Third Topo, shot ' num2str(shot)]);
        xlabel('Pixel Number (real space)');
        ylabel('Pixel Number (wavelength space)');
        set(h3, 'LineWidth', lnwdth);
        set(gca, 'FontSize', fntsz, 'XLim', [1 n_chan], 'YLim', [1 n_pix]);
        
        if saveTopo3 && savePlot
            fig_save = getframe(h1);
            [Xfig, mapfig] = frame2im(fig_save);
            imwrite(Xfig, [fileTopo3 '.png']);
        end
    end
    
    %% Plot Topo4
    if plotTopo4
        h1 = figure('Visible','on','Name','Topo 4','Position',...
            [S(3)/4, S(4)/6, S(3)/2 2*S(4)/3], 'Color', [1 1 1]);
        h2 = axes('Parent', h1, 'Position', [.1 .1 .8 .8], 'FontSize', fntsz);
        h3 = surf(X, Y, Ak(4)*squeeze(topos(4, :, :)));
        hold on;
        shading interp;
        colormap jet;
        colorbar;
        grid on;
        view([0 90]);
        title(['Fourth Topo, shot ' num2str(shot)]);
        xlabel('Pixel Number (real space)');
        ylabel('Pixel Number (wavelength space)');
        set(h3, 'LineWidth', lnwdth);
        set(gca, 'FontSize', fntsz, 'XLim', [1 n_chan], 'YLim', [1 n_pix]);
        
        if saveTopo4 && savePlot
            fig_save = getframe(h1);
            [Xfig, mapfig] = frame2im(fig_save);
            imwrite(Xfig, [fileTopo4 '.png']);
        end
    end
    
    %% Plot Topo5
    if plotTopo5
        h1 = figure('Visible','on','Name','Topo 5','Position',...
            [S(3)/4, S(4)/6, S(3)/2 2*S(4)/3], 'Color', [1 1 1]);
        h2 = axes('Parent', h1, 'Position', [.1 .1 .8 .8], 'FontSize', fntsz);
        h3 = surf(X, Y, Ak(5)*squeeze(topos(5, :, :)));
        hold on;
        shading interp;
        colormap jet;
        colorbar;
        grid on;
        view([0 90]);
        title(['Fifth Topo, shot ' num2str(shot)]);
        xlabel('Pixel Number (real space)');
        ylabel('Pixel Number (wavelength space)');
        set(h3, 'LineWidth', lnwdth);
        set(gca, 'FontSize', fntsz, 'XLim', [1 n_chan], 'YLim', [1 n_pix]);
        
        if saveTopo5 && savePlot
            fig_save = getframe(h1);
            [Xfig, mapfig] = frame2im(fig_save);
            imwrite(Xfig, [fileTopo5 '.png']);
        end
    end
    
    %% Plot Topo Combination
    if plotTopoComb
        h1 = figure('Visible','on','Name','First Three Topos','Position',...
            [S(3)/4, S(4)/6, S(3)/2 2*S(4)/3], 'Color', [1 1 1]);
        
        h2 = axes('Parent', h1, 'Position', [.12 .65 .8 .25], 'FontSize', fntsz);
        h3 = surf(X, Y, Ak(1)*squeeze(topos(1, :, :)));
        hold on;
        shading interp;
        colormap jet;
        colorbar;
        grid on;
        view([0 90]);
        title(['First 3 Topos, shot ' num2str(shot)]);
        set(h3, 'LineWidth', lnwdth);
        set(gca, 'FontSize', fntsz, 'XLim', [1 n_chan], 'YLim', combPixLim);
        set(gca, 'XTick', []);
        ylabel('A_1 * Topo 1');
        
        h2 = axes('Parent', h1, 'Position', [.12 .37 .8 .25], 'FontSize', fntsz);
        h3 = surf(X, Y, Ak(2)*squeeze(topos(2, :, :)));
        hold on;
        shading interp;
        colormap jet;
        colorbar;
        grid on;
        view([0 90]);
        set(h3, 'LineWidth', lnwdth);
        set(gca, 'FontSize', fntsz, 'XLim', [1 n_chan], 'YLim', combPixLim);
        set(gca, 'XTick', []);
        ylabel('A_2 * Topo 2');
        
        h2 = axes('Parent', h1, 'Position', [.12 .1 .8 .25], 'FontSize', fntsz);
        h3 = surf(X, Y, Ak(3)*squeeze(topos(3, :, :)));
        hold on;
        shading interp;
        colormap jet;
        colorbar;
        grid on;
        view([0 90]);
        set(h3, 'LineWidth', lnwdth);
        set(gca, 'FontSize', fntsz, 'XLim', [1 n_chan], 'YLim', combPixLim);
        ylabel('A_3 * Topo 3');
        xlabel('Pixel Number (real space)');
        
        if saveTopoComb && savePlot
            fig_save = getframe(h1);
            [Xfig, mapfig] = frame2im(fig_save);
            imwrite(Xfig, [fileTopoComb '.png']);
        end
    end
    
    %% Plot Chronos
    if plotChrons
        h1 = figure('Visible','on','Name','Topo 1','Position',...
            [S(3)/4, S(4)/6, S(3)/2 2*S(4)/3], 'Color', [1 1 1]);
        h2 = axes('Parent', h1, 'FontSize', fntsz);
        h3 = plot(time, V(:, 1), '-b');
        hold on;
        set(h3, 'LineWidth', lnwdth);
        h3 = plot(time, V(:, 2), '-r');
        set(h3, 'LineWidth', lnwdth);
        h3 = plot(time, V(:, 3), '-g');
        set(h3, 'LineWidth', lnwdth);
        h3 = plot(time, V(:, 4), '-c');
        set(h3, 'LineWidth', lnwdth);
        h3 = plot(time, V(:, 5), '-m');
        set(h3, 'LineWidth', lnwdth);
        
        title(['Chronos, shot ' num2str(shot)]);
        xlabel('time [ms]');
        ylabel('Normalized Chrono');
        set(h3, 'LineWidth', lnwdth);
        %     set(gca, 'FontSize', fntsz, 'XLim', [1 n_chan], 'YLim', [1 n_pix]);
        legend('1', '2', '3', '4', '5');
        
        if saveChrons && savePlot
            fig_save = getframe(h1);
            [Xfig, mapfig] = frame2im(fig_save);
            imwrite(Xfig, [fileChrons '.png']);
        end
    end
    
    
    %% Plot Chronos Weighted
    if plotChronsW
        h1 = figure('Visible','on','Name','Topo 1','Position',...
            [S(3)/4, S(4)/6, S(3)/2 2*S(4)/3], 'Color', [1 1 1]);
        h2 = axes('Parent', h1, 'FontSize', fntsz);
        h3 = plot(time, Ak(1)*V(:, 1), '-b');
        hold on;
        set(h3, 'LineWidth', lnwdth);
        h3 = plot(time, Ak(2)*V(:, 2), '-r');
        set(h3, 'LineWidth', lnwdth);
        h3 = plot(time, Ak(3)*V(:, 3), '-g');
        set(h3, 'LineWidth', lnwdth);
        h3 = plot(time, Ak(4)*V(:, 4), '-c');
        set(h3, 'LineWidth', lnwdth);
        
        title(['Weighted Chronos, shot ' num2str(shot)]);
        xlabel('time [ms]');
        ylabel('A_k * Chrono');
        set(h3, 'LineWidth', lnwdth);
        set(gca, 'FontSize', fntsz);
        legend('1', '2', '3', '4', 'Location', 'East');
        
        if saveChronsW && savePlot
            fig_save = getframe(h1);
            [Xfig, mapfig] = frame2im(fig_save);
            imwrite(Xfig, [fileChronsW '.png']);
        end
    end
    
    %% Save Data for Recombining
    if m == 1
        Ak_t = Ak(modes); % save weights
        V_t = V(:, modes); % save chronos
        topos_t = topos(modes, :, :); % save topos
    else
        Ak_p = Ak(modes); % save weights
        V_p = V(:, modes); % save chronos
        topos_p = topos(modes, :, :); % save topos
    end
    clear Ak V topos
    
end
%% Recombine Selected Modes into False Shot

if recombine
    data = zeros(n_time_o, n_pix_o, n_chan_o);
    data2 = data; % initialize both to zeros
    for n = 1:length(modes)
        for m = 1:size(V_t, 1) % loop over time
            data2(m, lam_bound, chan_bound_t) = Ak_t(n) * V_t(m, n) * topos_t(n, :, :);
            data2(m, lam_bound, chan_bound_p) = Ak_p(n) * V_p(m, n) * topos_p(n, :, :);
        end
        data = data + data2;
    end
    clear data2 Ak_t V_t topos_t Ak_p V_tp topos_p;
    
    % Reflip Data so everything in consistent
    
    data = data(:, end:-1:1, end:-1:1);
    
    % Save
    
    cd('T:\PhantomMovies');
    save(['shot' falseShot '.mat'], 'data');
    save(['t' falseShot '.mat'], 'time');
    
    disp(['Saved to shot ' falseShot ', modes: ' num2str(modes)]);
end




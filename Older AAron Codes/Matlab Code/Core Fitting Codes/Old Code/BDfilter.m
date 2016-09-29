function [data, time] = BDfilter(shot, nBDmodes, timeBound, param)
% Aaron Hossack
%
% The '3' version does a better job- now, tight bounds are set around the
% toroidal and poloidal data, they are individually BD'd, and then placed
% onto an array of zeros.  Thus, the movie is the same size as the original
% so the channel numbers match.

plotFits = 1;
BDwing = 6; % number of points to either side of center in wavelength space to include in BD domain

modes = 1:nBDmodes; % array of mode numbers to recombine and add together

first = param.peaks(1, 2);
last = param.peaks(find(param.peaks(:, 1) <= 36, 1, 'last'), 2);
chan_bound_t = floor(first):ceil(last);

first = param.peaks(find(param.peaks(:, 1) >= 37, 1), 2);
last = param.peaks(end, 2);
chan_bound_p = floor(first):ceil(last);

lam_bound = round(param.Center) - BDwing:round(param.Center) + BDwing;

%%
% Load Data Array

cd('T:\PhantomMovies');
data = importdata(['shot' int2str(shot) '.mat']); % [counts] (time x wavelength space x channel space)
time = importdata(['t' int2str(shot) '.mat']); % [ms]
cd('S:\MC_IDS\Matlab Code\Core Codes');

% Unflip images to correct for flip during Python conversion

data = data(:, end:-1:1, end:-1:1);

% BD specific time bound
if isempty(timeBound)
    timeBound = 1:length(time);
end
data = data(timeBound, :, :);
time = time(timeBound);

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
    else % poloidal fiber
        chan_bound = chan_bound_p;
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
    
    if plotFits
        disp(['Min = ' num2str(minVal)]);
        disp(['Median = ' num2str(median(BDdat(:)))]);
        disp(['Mean = ' num2str(mean(BDdat(:)))]);
        disp(['Std = ' num2str(std(BDdat(:)))]);
    end
    
    BDdat = BDdat - minVal; % subtract minimum value in data
    
    %% Perform SVD -----------------------
    
    [U, S, V] = svd(BDdat, 'econ');
    
    Ak = diag(S);
    
    % Rearrange data as it was before
    
    topos = zeros(n_time, n_pix, n_chan);
    
    for n = 1:n_time
        topos(n, :, :) = reshape(U(:, n), 1, n_pix, n_chan);
    end
    
    %% Plotting
    if plotFits
        [X, Y] = meshgrid(1:n_chan, 1:n_pix);
        
        S1 = get(0,'ScreenSize');
        lnwdth = 2;
        fntsz = 20;
        
        h1 = figure('Visible','on','Name','Topo 1','Position',...
            [S1(3)/4, S1(4)/6, S1(3)/2 2*S1(4)/3], 'Color', [1 1 1]);
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
    
end




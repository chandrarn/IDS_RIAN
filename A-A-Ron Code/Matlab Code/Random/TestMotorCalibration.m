% Calibration script for Hanao's fiber optics on HIT-SI

clear all; close all; clc;
%% Settings
testing = 0;
fileName = 'cal140825'; % file name that this script creates
motorShot = 14082501; % shot with motor calibration
bins = [38, 51;...
    69, 87;...
    103, 122;...
    139, 155;...
    172, 190;...
    206, 225;...
    240, 258;...
    276, 292];

%shots = [14082901:14082907]; % stationary shots, fibers 1 & 2, 1 & 3, etc..., 1 & 8
shots = [14082503:14082509]; % stationary shots, fibers 1 & 2, 1 & 3, etc..., 1 & 8

motorSpeed = 1; % [nm/s]

%% Load Data
addpath('T:\PhantomMovies');
data = importdata(['shot' int2str(motorShot) '.mat']); % [counts] (time x wavelength space x channel space)
time = importdata(['t' int2str(motorShot) '.mat']); % [ms]

data = cast(data, 'double');
data = data(:, end:-1:1, end:-1:1);

[n_time, n_wavelength, n_spatial] = size(data);

% Bin Data

for n = 1:n_time
    for m = 1:size(bins, 1)
        dataB(n, m, :) = sum(squeeze(data(n, :, bins(m, 1):bins(m, 2))), 2);
    end
%     if testing
%         plot(squeeze(dataB(n, 1, :)));
%         hold all;
%         for m = 2:size(bins, 1)
%             plot(squeeze(dataB(n, m, :)));
%         end
%         pause(2);
%         clf;
%     end
end

% Data is now binned [n_time, n_chan (8), n_wavelength (96)]

%% Fit Gaussians from Motor Shot
options = optimsetv61('lsqcurvefit'); % set options for curve fitting
options.TolFun = 1e-8 * max(dataB(:)); % set tolerance for curve fitting
f = @singletFun; % function handle

motorCenters = zeros(n_time, size(bins, 1));
xdata = 1:n_wavelength; % 1:96
for m = 1:size(bins, 1) % Loop over channels
    for n = 1:n_time % loop over time
        ydata = squeeze(dataB(n, m, :))'; % curve for one channel, one time
        
        [amp, x0] = max(ydata);
        fwhm = length(find(ydata > 0.5 * amp));
        sigma = fwhm / 2.35;
        area = 0.8 * sqrt(2*pi) * sigma * amp; % 0.8 factor arbitrary for best initial guess
        offset = min(ydata);
        
        guess = [area, x0, sigma, offset];
        lb = [0.3*area, guess(2)-1, 0.2*guess(3), 0]; % lower bound
        ub = [3*area, guess(2)+1, 3*guess(3), 0.5*amp]; % upper bound
        parameters = lsqcurvefit(f, guess, xdata, ydata, lb, ub, options);
        
        motorCenters(n, m) = parameters(2);
%         if testing
%             xfine = linspace(xdata(1), xdata(end), 1000);
%             yfine = singletFun(parameters, xfine);
%             yguess = singletFun(guess, xfine);
%             
%             plot(xfine, yfine, '-b'); % fit evaluated
%             hold on;
%             plot(xfine, yguess, '-c'); % initial guess
%             plot(squeeze(dataB(n, m, :)), '+r'); % raw data points
%             plot([motorCenters(n, m), motorCenters(n, m)], [0, amp], '-k'); % line marking center
%             title(['Channel ' num2str(m) ', Time ' num2str(n)]);
%             xlabel('Wavelength Direction [pixels]');
%             ylabel('Amplitude [counts]');
%             pause(1);
%             clf;
%         end
    end
end
    
%% Calculate PIX_SP
nm = motorSpeed * 1e-3 * time; % convert 'time' to [s], [s] to [nm]

for m = 1:size(bins, 1)
    p = polyfit(motorCenters(:, m), nm, 1);
    PIX_SP(m) = abs(p(1));
end
PIX_SP = 1e-9 * PIX_SP; % convert from nm to meters
disp(num2str(PIX_SP));

clear dataB data

%% Load Data for REL_INT and SIGMA

for m = 1:length(shots) % 1:7
    data = importdata(['shot' int2str(shots(m)) '.mat']); % [counts] (time x wavelength space x channel space)
%     time = importdata(['t' int2str(motorShot) '.mat']); % [ms]
    
    data = cast(data, 'double');
    data = data(:, end:-1:1, end:-1:1);
    
    [n_time, n_wavelength, n_spatial] = size(data);
    
    % Time Average Data
    data = squeeze(mean(data, 1));
    
    %% Channel 1
    ydata = sum(data(:, bins(1, 1):bins(1, 2)), 2)';
    
    [amp, x0] = max(ydata);
    fwhm = length(find(ydata > 0.5 * amp));
    sigma = fwhm / 2.35;
    area = 0.8 * sqrt(2*pi) * sigma * amp; % 0.8 factor arbitrary for best initial guess
    offset = min(ydata);
    
    guess = [area, x0, sigma, offset];
    lb = [0.3*area, guess(2)-1, 0.2*guess(3), 0]; % lower bound
    ub = [3*area, guess(2)+1, 3*guess(3), 0.5*amp]; % upper bound
    parameters = lsqcurvefit(f, guess, xdata, ydata, lb, ub, options);
    ch1_int = parameters(1);
    
    if m == 1
        SIGMA(m) = parameters(3);
        REL_INT(m) = 1;
        CENTER(m) = parameters(2);
    end
    
    if testing
        xfine = linspace(xdata(1), xdata(end), 1000);
        yfine = singletFun(parameters, xfine);
        yguess = singletFun(guess, xfine);
        
        plot(xfine, yfine, '-b'); % fit evaluated
        hold on;
        plot(xfine, yguess, '-c'); % initial guess
        plot(xdata, ydata, '+r'); % raw data points
        title(['Channel ' num2str(m)]);
        xlabel('Wavelength Direction [pixels]');
        ylabel('Amplitude [counts]');
        pause(1);
        clf;
    end
    
    %% All Other Channels
    
    ydata = sum(data(:, bins(m+1, 1):bins(m+1, 2)), 2)'; % Bin data
    
    [amp, x0] = max(ydata);
    fwhm = length(find(ydata > 0.5 * amp));
    sigma = fwhm / 2.35;
    area = 0.8 * sqrt(2*pi) * sigma * amp; % 0.8 factor arbitrary for best initial guess
    offset = min(ydata);
    
    guess = [area, x0, sigma, offset];
    lb = [0.3*area, guess(2)-1, 0.2*guess(3), 0]; % lower bound
    ub = [3*area, guess(2)+1, 3*guess(3), 0.5*amp]; % upper bound
    parameters = lsqcurvefit(f, guess, xdata, ydata, lb, ub, options);
    
    SIGMA(m+1) = parameters(3);
    REL_INT(m+1) = parameters(1) / ch1_int;
    CENTER(m+1) = parameters(2);
    
end

disp(num2str(SIGMA));
disp(num2str(REL_INT));
    
if ~testing
    save(['T:\Hanao\Calibration\' fileName], 'PIX_SP', 'SIGMA', 'REL_INT', 'CENTER', 'bins');
end










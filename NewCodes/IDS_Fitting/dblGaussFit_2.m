function[fit_par, dPar, stddev, param] = dblGaussFit_2(data2, param, options)

% Revised 2-17-13 by ACH to implement the Levenberg-Marquart method of
% curve fitting and error estimation.
%
% Revised 2-21-13 by ACH to first curve fit with lsqcurvefit, then put
% final parameters into LM method to find errors.

addpath('S:\General Matlab\LM Method');

%% Initial parameter estimates

[n_time, ymax, xmax] = size(data2); % determine size of data array
                                    % [time, pixel, channel]

sumOfFrame = sum(sum(data2, 3), 2); % sums over two data dimensions of 'data2' array leaving
                                   % scalar sum of all counts for each time.
                                   
[dummy, n_bright] = max(sumOfFrame); % identify time index with most light

[dummy, chan_bright] = max(sum(data2(n_bright, :, :), 2)); % identify brightest channel

% find approximate offset

offset = min(data2(n_bright, :, chan_bright));

% Convert wavelength separation of doublet into channel separation

param.deltaPix = (param.LineLam(1) - param.LineLam(2)) ./ param.PIX_SP; % [m / (m/pix) = pix]

% Find approximate ratio of line intensities

oneDSpectr = squeeze(sum(sum(data2, 3), 1)); % collapse data in time and channel dimensions
                                             % to get 1-D spectrum
                                             % 'squeezed' into one column vector
                                             
                                             plot(oneDSpectr);
                                             
spectrOffst = min(oneDSpectr); % total 1-D spectrum offset
                                
[bright, indBright] = max(oneDSpectr); % find maximum intensity and pixel number (index)

dPix = median(param.deltaPix); % median pixel spacing for use in next operations

try
    dim1 = max(oneDSpectr(indBright + round(dPix) - 1 : indBright + round(dPix) + 1));
    % if the other line is a longer wavelength, look within a small range at
    % the expected separation distance
catch
    dim1 = 0; % if the code goes here the attempt to find the dim line at
              % longer wavelength failed because it looked out of bounds
end

try
[dim2, indDim2] = max(oneDSpectr(indBright - round(dPix) - 1 : indBright - round(dPix) + 1));
    % if the other line is a shorter wavelength, look within a small range at
    % the expected separation distance
indDim2 = indDim2 + indBright - round(dPix) - 2;

catch
    dim2 = 0; % if the code goes here the attempt to find the dim line at
              % shorter wavelength failed because it looked out of bounds
end

if dim1 > dim2 % the dim line is at a longer wavelength
    center = indBright; % the 'center' parameter used for velocity fitting
                        % is arbitrarily chosen in 'loadParams.m' to be the
                        % lower wavelength line, ergo for this situation
                        % the bright line is the lower wavelength one and
                        % is used for the center estimate
    amplitude = data2(n_bright, indBright, chan_bright); % first guess at primary line amplitude
    param.ratio = (dim1 - spectrOffst) / (bright - spectrOffst); % line amplitude ratio used for curve fitting
    
else % the dim line is at a shorter wavelength
    center = indBright - dPix; % initial guess for 'center' parameter, always the shorter wavelength
    amplitude = data2(n_bright, indDim2, chan_bright); % initial guess for primary line amplitude (always shorter wavelength)
    param.ratio = (bright - spectrOffst) / (dim2 - spectrOffst); % line amplitude ratio for curve fitting

end

% find approximate width

fwhm = 2 * sum(data2(n_bright, 1:round(center), chan_bright) - offset >= 0.5* (amplitude - offset)); % find crude FWHM
    % twice the number of pixel on the low wavelength side of the low
    % wavelength line with amplitude >= 1/2 the primary line

width = fwhm / 2.35482; % convert to standard Gaussian function parameter

area = amplitude * width * sqrt(2*pi); % using area of gaussian

guess = [area, center, width, offset]; % initial guess

%% Fit to Gaussians

% Manually set weight array for 'lm_sigma' error analysis
% weight = 7*ones(1, ymax)';
weight = 0;

% Prepare empty arrays
fit_par = NaN*zeros(n_time, 4, param.n_chan); % preallocate fit parameters
blank1 = NaN*ones(1, 4); % overwrite parameters if found to be bad
blank2 = NaN*ones(1, ymax); % overwrite stddev data if bad
stddev = NaN*zeros(n_time, ymax, xmax); % preallocate standard deviation
                                        % array same size as data array
dPar = fit_par; % uncertainty in fit parameters
dp = [0.001, 0.001, 0.001, 0.001]; % fractional increment of 'p' for numerical derivatives

xdata = 1:ymax; % vector of pixel numbers in wavelength space
% dxData = []; % These could be vectors of known uncertainties in x (pixels)
% dyData = []; %     "                  "               "         y (data values)
f2 = @doubletFunLM; % LM code must have slightly different form
                    % of function

for m = 1:param.n_chan % order of 'for' loops switched from normal convention
                       % so that 'f' is not redefined at every iteration
    c = [param.deltaPix(m), param.ratio]; % extra fitting constants
           
    f1 = @(guess, xdata)doubletFun(guess, xdata, param.deltaPix(m), param.ratio);
        % this is an "anonymous function", the purpose of which is to allow
        % the specification of extra unchanging parameters.
    
    for n = 1:n_time
        % Curve Fit using 'lsqcurvefit' first
        [fit_par(n, :, m), resnorm, residual, exitflag] = ...
            lsqcurvefit(f1, guess, xdata, data2(n, :, m), param.lb, param.ub, options);
        
        % exclude bad data
        if or(exitflag == 0, sqrt(resnorm) > fit_par(n, 1, m)*param.Resid_to_Sig)
            fit_par(n, :, m) = blank1;
            dPar(n, :, m) = blank1;
            stddev(n, :, m) = blank2;
        else
            % Find error estimates using modified LM code
            [X2, dPar(n, :, m), stddev(n, :, m), corr, R_sq] = ...
                lm_sigma(f2, fit_par(n, :, m), xdata, data2(n, :, m)', dp, c, weight);
        end


        % exclude bad data ---------  REVISIT THIS
%         if sqrt(resnorm(n, m)) > fit_par(n, 1, m)*param.Resid_to_Sig
%             fit_par(n, :, m) = blank;
%         end
        
    end
end
end
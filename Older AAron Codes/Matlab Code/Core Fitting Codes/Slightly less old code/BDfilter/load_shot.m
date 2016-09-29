function[out, param] = load_shot(shot, nBDmodes, timeBound, line, frame_sum, plots, s,hitsi3)
% This code handles all the analysis of a MultiChord IDS shot.  The only
% inputs it needs are the shot number, number of frames to sum
% (typically 1), time override, and peak locations.  If 't_override' is set to 0, the time
% to analyze will be determined by 'shotDuration.m', or 't_override' can be a
% two value vector with start and end times in ms.
addpath('C:\Program Files\MDSplus\MATLAB');
addpath('T:\RChandra\A-A-Ron Code\Matlab Code\Core Fitting Codes\lsqcurvefit');
cd('T:\RChandra\A-A-Ron Code\Matlab Code\Core Fitting Codes\BDfilter');

% Load Constants and Parameters

[param, options] = loadParams(shot, line,hitsi3);

%% Special Preparations

% BD filtering

if isempty(nBDmodes)
    % Load Data Array
    %Try to import the files from the new format, if possible
    try
        importdata(['Shot ' int2str(shot) '.mat']);
        data=CineArray;
        time=TimeVector;
        newCineType=1;
    catch
        data = importdata(['shot' int2str(shot) '.mat']); % [counts] (time x wavelength space x channel space)
        time = importdata(['t' int2str(shot) '.mat']); % [ms]
        newCineType=0
        
        % Unflip images to correct for flip during Python conversion
        data = data(:, end:-1:1, end:-1:1);
        data = cast(data, 'double'); % comes in as 'uint16' and ruins everything
        %Olny if you're NSTX and dont know whats up.
    end
 
else
    [data, time] = BDfilter(shot, nBDmodes, timeBound, param);
end
    

% Frame averaging, if applicable

% Note from June 2013 - these are probably all broken due to switch to 2D
% fitting
% Note from April 2014 - Thank you, A-A-Ron

if frame_sum ~= 1
    [time, data] = frameSum(time, data, frame_sum, param);
end

% Summing chords together

if s.spSum
    [data, param] = chordSum(data, param, s);
end

% Adding a pattern of time points together

if s.tSum
    [data, time, param] = timeSum(data, param, s);
end

%% Calculations

% Gaussian Fitting
tic
if param.dbl % fitting to a doublet
    [fit_par, dPar, bounds, stddev, param, guesses] = dblGaussFit_2D(data, param, options);
else % single line
    [fit_par, dPar, bounds, stddev, param, guesses] = gaussFit_2D(data, param, options);
end
toc
% Calculate Temperature, Velocity, Area

[temp, vel, int, param] = calcPhysics(fit_par, param);

%% Error
if param.calcError
    % Caulculate Upper Error Physics

    [tempU, velU, dummy, dummy] = calcPhysics(fit_par + dPar, param);

    % Caulculate Lower Error Physics

    [tempL, velL, dummy, dummy] = calcPhysics(fit_par - dPar, param);

    % Convert to Error Bars

    dtempU = tempU - temp;
    dvelU = velU - vel;
    dtempL = temp - tempL;
    dvelL = vel - velL;
end

%% Output

% Sanitize output by discarding ridiculous data

[temp, vel] = sanitizeOut(temp, vel, param.limits);

% Bundle all outputs as structure

out.time = time; % (n_time)
out.temp = temp; % (n_time) x (param.n_chan)
out.vel = vel; % (n_time) x (param.n_chan)
out.int = int; % intensity, relatively calibrated (n_time) x (param.n_chan)
% out.resnorm = resnorm; % (n_time) x (param.n_chan)
out.fit_par = fit_par; % (n_time) x (param.n_chan) x (6 (number of parameters))
out.data = data; % (n_time) x (pixels in wavelength space) x (pixels in real space)
out.bounds = bounds; % (n_time) x (param.n_chan) x (4: x1, x2, y1, y2)
out.guesses = guesses; % (n_time) x (param.n_chan) x (6 (number of parameters))

if param.calcError
    % Added when switched to Levenberg-Marquardt error method
    out.dparams = dPar; % same size as 'out.params'- standard error in params
    out.stddev = stddev; % standard deviation for every data point
    out.tempU = dtempU; % upper error bar for temperature
    out.tempL = dtempL; % lower error bar for temperature
    out.velU = dvelU; % upper error bar for velocity
    out.velL = dvelL; % lower error bar for velocity
end

end



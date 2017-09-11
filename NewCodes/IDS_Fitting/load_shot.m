function[out, param] = load_shot(shot, nBDmodes, timeBound, line, frame_sum, plots, s, hitsi3, useTree)
% This code handles all the analysis of a MultiChord IDS shot.  The only
% inputs it needs are the shot number, number of frames to sum
% (typically 1), time override, and peak locations.  If 't_override' is set to 0, the time
% to analyze will be determined by 'shotDuration.m', or 't_override' can be a
% two value vector with start and end times in ms.
addpath('C:\Program Files\MDSplus\MATLAB');
%addpath('T:\RChandra\A-A-Ron Code\Matlab Code\Core Fitting Codes\lsqcurvefit');
%cd('T:\RChandra\A-A-Ron Code\Matlab Code\Core Fitting Codes\BDfilter');
addpath('T:\PhantomMovies\');


% Load Constants and Parameters

[param, options] = loadParams(shot, line, hitsi3, useTree);
assignin('base','param',param)

%% Special Preparations

size(param.peaks)
%% BD filtering

if isempty(nBDmodes)
    % Load Data Array
    %Try to import the files from the new format, if possible
        try
            dummy = importdata(['Shot ' int2str(shot) '.mat']);
            data = dummy.CineArray;
            time = dummy.TimeVector;
            % Sometimes PhantomStalker appears to flip the y-axis, should be fixed
            %data = data(end:-1:1, :,:);%precut time
            data=shiftdim(data,2); % Necessary to match NSTX codes
            newCineType = 1;
            clear dummy;
        catch
            if ~isempty(s.sim)
                shot = str2num([int2str(s.sim) int2str(shot)]); % FUCK YEAH
            end
            data = importdata(['shot' int2str(shot) '.mat']); % [counts] (time x wavelength space x channel space)
            time = importdata(['t' int2str(shot) '.mat']); % [ms]
            newCineType = 0;

            % Unflip images to correct for flip during Python conversion
            data = data(:, end:-1:1, end:-1:1);
            data = cast(data, 'double'); % comes in as 'uint16' and ruins everything
            %Olny if you're NSTX and dont know whats up.
        end
        data = data(timeBound,:,:);
        time = time(timeBound);
        figure; 
        surf(squeeze(sum(data, 1)./size(data,1))); 
        title(['Shot: ' num2str(shot)]);
        ylabel('\lambda Axis');
        xlabel('Spatial Axis');
        shading interp;
        view([ 0 -90]);
    
else
    [data, time] = BDfilter(shot, nBDmodes, timeBound, param, s);
end

assignin('base','data',data);
assignin('base','options',options);

%% Frame averaging, if applicable

% Note from June 2013 - these are probably all broken due to switch to 2D
% fitting
% Note from April 2014 - Thank you, A-A-Ron

if frame_sum ~= 1
    [time, data] = frameSum(time, data, frame_sum, param);
end

%% Summing chords together

if s.spSum
    [data, param] = chordSum(data, param, s);
end

%% Adding a pattern of time points together

if s.tSum
    [data, time, param] = timeSum(data, param, s);
end

%% Calculations

% Gaussian Fitting

%% Temp for fixing old data
%dat = importdata(['T:\IDS\Data Repository\dat' num2str(shot) '10.mat']);

for nn = 1:length(line)
    tic
    [fit_par, dPar, bounds, stddev, param, guesses] = gaussFit_2D(data, param, options, nn);
    %fit_par = dat(nn).fit_par;
    %bounds = dat(nn).bounds;
    %guesses = dat(nn).guesses;
    %stddev = NaN;
    %dPar = NaN;
    
    
    toc
    
    % Calculate Temperature, Velocity, Area
    [temp, vel, int, param] = calcPhysics(fit_par, param, s, nn);
    
    %% Error
    if param.calcError
        % Caulculate Upper Error Physics
        
        [tempU, velU, intU, ~] = calcPhysics(fit_par + dPar, param, s, nn);
        
        % Caulculate Lower Error Physics
        
        [tempL, velL, intL, ~] = calcPhysics(fit_par - dPar, param, s, nn);
        
        % Convert to Error Bars
        
        dintU = intU - int;
        dtempU = tempU - temp;
        dvelU = velU - vel;
        dintL = int - intL;
        dtempL = temp - tempL;
        dvelL = vel - velL;
    end
    
    %% Output
    
    % Sanitize output by discarding ridiculous data
    assignin('base','PreSanitizedDataTemp',temp);
    assignin('base','PreSanitizedDataVel',vel);
    [temp, vel] = sanitizeOut(temp, vel, param.limits);
    
    % Bundle all outputs as structure
    
    out(nn).time = time; % (n_time)
    out(nn).temp = temp; % (n_time) x (param.n_chan)
    out(nn).vel = vel; % (n_time) x (param.n_chan)
    out(nn).int = int; % intensity, relatively calibrated (n_time) x (param.n_chan)
    % out.resnorm = resnorm; % (n_time) x (param.n_chan)
    out(nn).fit_par = fit_par; % (n_time) x (param.n_chan) x (6 (number of parameters))
    out(nn).data = data; % (n_time) x (pixels in wavelength space) x (pixels in real space)
    out(nn).bounds = bounds; % (n_time) x (param.n_chan) x (4: x1, x2, y1, y2)
    out(nn).guesses = guesses; % (n_time) x (param.n_chan) x (6 (number of parameters))
    
    if param.calcError
        % Added when switched to Levenberg-Marquardt error method
        out(nn).dparams = dPar; % same size as 'out.params'- standard error in params
        out(nn).stddev = stddev; % standard deviation for every data point
        out(nn).intU = dintU; % upper error bar for intensity
        out(nn).intL = dintL; % lower error bar for intensity
        out(nn).tempU = dtempU; % upper error bar for temperature
        out(nn).tempL = dtempL; % lower error bar for temperature
        out(nn).velU = dvelU; % upper error bar for velocity
        out(nn).velL = dvelL; % lower error bar for velocity
    end
end
assignin('base','out',out);
end



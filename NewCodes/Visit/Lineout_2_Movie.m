%% Lineout_2_Movie.m
%
% This script takes the lineout files generated in VisIt/Python in the
% format 'LO <tet/nim> <tor=1/pol=0> _ <config> _ <time point> .mat', ie:
% LO1_3_63.mat and generates a false movie as if it had been created by the
% fast camera.
%
% The script loops through all times and both fibers, if applicable.
%
% 'lineVals{p, 1}.vx' contains the pth channel and velocity in the x
% direction. Channel numbers are matched to expermental channel numbers in
% 'chan_range'
%
% CONFIG REFERENCE:
% 1 toroidal, 71 degree port
% 2 toroidal, mohawk port in midplane
% 3 toroidal, axial port at 135 degrees
% 4 toroidal, mohawk port perp.
% 1 axial fiber @ toroidal angle 45
% 2 axial fiber @ toroidal angle 270
% 3 axial fiber @ toroidal angle 135
clear all; close all; clc;
addpath('~/IDS/Matlab/');
addAllThePaths;

%% INPUT SETTINGS
testing = 0; % make lots of plots along the way
saveTemp = 0; % special setting for outputting weighted, chord averaged temperature.

refShot = 12949910; % real IDS shot including BD number added on
refLine = 3; % corresponds to 'line' setting in Batch_Correct and the tree, typically 3 for C III triplet
refLineDat = 2; % the index of the above line, ie: dat(#).temp, typically 2 for C III triplet
sim = 'tet'; % 'tet' or 'nim'
nimDataSet = 'fin_beta_ramp'; % folder containing everything related to a particular run
nimMatFile = '24x24_14_injRamp'; % for creating time base
tetDataSet = 'aaronData_141022'; % folder containing everything related to a particular run
tetMatFile = 'PSITET_14kHz_xMHD'; % for creating time base
tor = 1; % tor = 1, pol = 0 (refers to fiber)
config = 2; % config of the primary dataset

tor2 = []; % leave blank '[]' if only using one fiber
config2 = 1; % config of the secondary dataset

newShotNum = '5129499'; % new shot number, leading number indicates type of simulation, followed by real shot number
                      % 2 = NIMROD --- weighted by n^2 ---  NO temperature  --- 7x timeSum
                      % 3 = PSI-TET --- weighted by n^2 ---  NO temperature  --- 7x timeSum
                      % 4 = NIMROD --- weighted by n^2 ---  NO temperature  --- no timeSum
                      % 5 = PSI-TET --- weighted by n^2 --- with temperature --- 7x timeSum
                      % 6 = NIMROD --- weighted by n^2 --- with temperature --- 7x timeSum
                      % 7 = PSI-TET -- weighted by n^2 --- with temperature --- no timeSum
                      % 8 = NIMROD --- weighted by n^2 --- with temperature --- no timeSum

timeSum = [7]; % number of points to add together
              % leave blank if none
                      
chWeight = 1; % 0 = no chord weighting
              % 1 = n^2

tempBroaden = 1; % 0 = do not use temperature
                 % 1 = broaden lines as a function of temperature

denFac = 1e-19; % factor to multiply by density to keep numbers down
dataFac = 0.07; % scale factor for final data array to make values similar to IDS

%% LOAD IDS DATA, EXTRACT INFO

load(['dat' int2str(refShot) '.mat']);

% Extract settings from IDS 'mat' file
[~, yPix, xPix] = size(dat(1).raw); % size of original IDS movie

% Extract info for "binning" velocity data
% first, find velocity per pixel:
vPerPix = dat(1).param.c / dat(1).param.LineLam(refLine) * dat(1).param.PIX_SP;

%% LOAD LINEOUTS

if strcmp(sim, 'tet')
    simFolder = 'PsitetData';
    simDataSet = tetDataSet;
    simMatFile = tetMatFile;
elseif strcmp(sim, 'nim')
    simFolder = 'NimrodData';
    simDataSet = nimDataSet;
    simMatFile = nimMatFile;
end
cd(['/home/aaron/IDS/' simFolder '/' simDataSet '/LO/']);

matList = dir(['LO' sim int2str(tor) '_' int2str(config) '_*.mat']); % all times for this fiber and configuration

%% Calculate Normal Vectors Along Chords
load(matList(1, 1).name); % loads 'lineVals' and 'chan_range'
simDat.peaks = chan_range; % for feeding the channels to 'findIDSforVisIt' which it interprets as dat.peaks
ptShift = 0;
lengthBeyond = 0.1;
[~, vec1] = findIDSforVisIt(tor, config, simDat, ptShift, lengthBeyond); % first fiber

%% Second Fiber
if ~isempty(tor2)
    matList2 = dir(['LO' sim int2str(tor2) '_' int2str(config2) '_*.mat']); % all times for other fiber and configuration
    matList = horzcat(matList, matList2);
    load(matList(1, 2).name); % loads 'lineVals' and 'chan_range'
    simDat.peaks = chan_range; % for feeding the channels to 'findIDSforVisIt' which it interprets as dat.peaks
    [~, vec2] = findIDSforVisIt(tor2, config2, simDat, ptShift, lengthBeyond); % optional second fiber
end

%% PREPARE SETTINGS

nTimes = size(matList, 1);
data = zeros(nTimes, yPix, xPix); % initialize 'movie' data array

if saveTemp
    avgTemp = zeros(nTimes, length(chan_range));
end

for m = 1:nTimes
    for k = 1:size(matList, 2) % 1 or 2 if using both fibers
        load(matList(m, k).name); % loads 'lineVals' and 'chan_range'
        
        %% Calculate v along the chords
        if k == 1
            vec = vec1;
        else
            vec = vec2;
        end
        
        for p = 1:length(chan_range)
            % dot velocity vetor onto unit chord vector
            % !!! minus sign put in because of velocity convention. 'vec'
            % is a unit vector pointing toward the telescope, but the
            % prevailing theory is positive velocity points away.
            vdl(p).v = -[lineVals{1, p}.vx(2:2:end)', lineVals{1, p}.vy(2:2:end)', lineVals{1, p}.vz(2:2:end)'] * vec(p, :)'; 
        end
        %
        if testing
            figure(eval([num2str(m) num2str(k) '1']))
            hold all;
            for p = 1:length(chan_range)
                plot(vdl(p).v);
            end
            title('v \cdot dl');
            xlabel('Index');
            ylabel('[m/s]');
        end
        %        
        
        %% Chord Weighting
        switch chWeight
            case 0 % no weighting
                for p = 1:length(chan_range)
                    wVec(p).w = ones(size(vdl(p).v, 1), size(vdl(p).v, 2));
                end
                
            case 1 % n^2
                for p = 1:length(chan_range) % normalize density by a factor
                    wVec(p).w = (denFac * lineVals{1, p}.n(2:2:end)).^2; % (density, n, scaled) squared
                end
                
        end
        %{
        if testing
            figure(eval([num2str(m) num2str(k) '2']))
            hold all;
            for p = 1:length(chan_range)
                plot(wVec(p).w);
            end
            title('Chord Weight');
            xlabel('Index');
            ylabel('Fractional Weight');
        end
        %}
        
        for p = 1:length(chan_range)
            % find real channel index corresponding to 'chan_range(p)'
            dataInd = find(dat(1).peaks == chan_range(p));
            
            if ~isempty(dataInd) % I have calibration data for this channel
                %% Create 2D Grid Identical to Real Data
                x0 = round(dat(1).param.peaks(dataInd, 2));
                y0 = round(dat(1).param.Center(dataInd, refLineDat));
        
                xBound = x0 - dat(1).param.xWing : x0 + dat(1).param.xWing;
                yBound = y0 - dat(1).param.yWing : y0 + dat(1).param.yWing;
                
                [X, Y] = meshgrid(xBound, yBound);
                
                % Reshape Grid
                x(:, 1) = X(:);
                x(:, 2) = Y(:);
                
                % Evaluate the 2D temperature
                % Gaussian at every velocity point independently along
                % the chord. Sum all the Gaussians with optional
                % density weighting to get total distribution.
                %
                % NB: this covers temperature AND velocity, so only
                % need to convolve this with inst. function.
                
                if tempBroaden
                    root_T = sqrt(lineVals{1, p}.T(2:2:end)); % prepare sqrt(T)

                    % All non-changing constants to convert sqrt(T) into
                    % sig_y. sig_y = root_T * sig_coeff.
                    sig_coeff = sqrt(dat(1).param.LineLam(refLine)^2 * dat(1).param.kBoltz * 11605 / ...
                        (dat(1).param.c^2 * dat(1).param.IonMass(refLineDat) * dat(1).param.PIX_SP(dataInd)^2));
                    % par_5 = sqrt( sig_y(temp)^2 + sig_y(inst)^2 )
                    par_5 = sqrt((sig_coeff * root_T).^2 + dat(1).param.peaks(dataInd, 5).^2); % [vector]
                else
                    par_5 = dat(1).param.peaks(dataInd, 5) * ones(1, length(vdl(p).v)); % [vector]
                end
                
                % calculate "volume" for all points ahead of time
                % still need to do: par(1) = par_1(w) * sig_y
                par_1 = 2 * pi * dat(1).param.peaks(dataInd, 4) * wVec(p).w .* par_5; % [vector]
                
%                 par_2 = dat(1).param.Center(dataInd, refLineDat); % x_0, [scalar]
                par_2 = dat(1).param.peaks(dataInd, 2); % x_0, [scalar]
                
                % par 3 (y_0) is every velocity point along the chord.
                % this converts velocity to pixels, offsets so 0 matches experimental zero.
                par_3 = dat(1).param.Center(dataInd, refLineDat) + (1 / vPerPix(dataInd)) * vdl(p).v; % [vector]
                
                par_4 = dat(1).param.peaks(dataInd, 4); % sig_x [scalar]
                par_6 = 0; % no offset
                
                z = zeros(size(x, 1), 1); % initialize empty array
                for q = 1:length(vdl(p).v) % loop over all velocity points for this chord/channel
                    z = z + singletGauss2D([par_1(q), par_2, par_3(q), par_4, par_5(q), par_6], x); % evaluate function
                end
                Z = reshape(z, size(X, 1), size(X, 2));
                
                if testing
                    figure(eval([num2str(m) num2str(k) '8']))
                    clf;
                    surf(X, Y, Z);
                    hold on;
                    view([0 90]);
                    colorbar;
                    shading interp;
                    title('Convolved Inst. Func., Vel. Dist. and Temp.');
                end
                
                if saveTemp
                    avgTemp(m, p) = sum(wVec(p).w .* lineVals{1, p}.T(2:2:end)) / sum(wVec(p).w);
                end
                
            else
            end
            
            %% Add New Channel Data to 'data' array
            data(m, yBound, xBound) = Z;
            
        end % channel loop
    end % fiber loop
end % time loop

%% Create Time Base
cd(['/home/aaron/IDS/' simFolder '/' simDataSet '/DUMP/']);

t = []; % initialize
if strcmp(sim, 'tet')
    list = dir('out*.xmf');
    for n = 1:length(list)
        text = fileread(list(n).name);
        tmp = regexp(regexp(text,'Time Value="[0-9.]*"','match'),'[0-9.]*','match');
        time = str2num(char(tmp{1}));
        t = horzcat(t, time);
    end
elseif strcmp(sim, 'nim')
    load(['/home/aaron/IDS/' simFolder '/' simDataSet '/' simMatFile]);
    for n = 1:length(nimsave.times)
        if exist(['dump_' sprintf('%05d', nimsave.times(n)) '_b.vtk'], 'file')
            t = [t, nimsave.time(n)];
        end
    end
end

t = 1e3 * t; % convert from s to ms

%% Time Averaging

if ~isempty(timeSum)
    n_time2 = floor(size(data, 1) / timeSum);
    data2 = zeros(n_time2, size(data, 2), size(data, 3));
    t2 = zeros(1, n_time2);
    for n = 1:n_time2
        data2(n, :, :) = sum(data((n-1) * timeSum + 1 : n * timeSum, :, :), 1);
        t2(n) = mean(t((n-1) * timeSum + 1 : n * timeSum));
    end
    data = data2;
    t = t2;
    
    dataFac = dataFac / timeSum; % rescale so data is still a reasonable value
end

%% Save Movie

% Throw in scale factor to make lines more closely match HIT-SI data pixel
% values
data = dataFac * data;

% Reflip Data so everything in consistent
data = data(:, end:-1:1, end:-1:1);

if ~testing
    cd(['/home/aaron/IDS/' simFolder '/' simDataSet]);
    save(['shot' newShotNum '.mat'], 'data');
    save(['t' newShotNum '.mat'], 't');
    if saveTemp
        save(['temp' newShotNum '.mat'], 'avgTemp');
    end
end




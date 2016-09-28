% Aaron Hossack
% Dec. 4th, 2013
% 
% First attempt at putting IDS data into .vtk format for VisIt
%
clear all; close all; clc;
addpath('S:\MC_IDS\Matlab Code\Data Repository');
addpath('S:\General Matlab');
addpath('S:\MC_IDS\Matlab Code\Display'); % needed for 'trimRange'
addpath('S:\MC_IDS'); % where NIMROD poloidal boundary is
addpath('S:'); % where text files containing line out velocity are

%% Settings

shot = 12949910; % shot number corresponding to 'dat' file
% shiftVel = -7000; % [m/s] velocity shift
%                   % -7000 for 129499 and similar
%                   % 16000 for 129591 and similar
inDir = 'S:\'; % directory where input file is
inFile = 'Line_nim_dynamic_min';
tor = 1; % 1 = toroidal array (1:36), 0 = poloidal array (37:72)
outDir = 'S:\MC_IDS\Psi TET Comparisons\Output\';
str = date;
comment = ['BD, ' str];

configt = 2;
configp = 1;

%% Load and Prepare Data

load(['dat' num2str(shot)]);

data = importdata([inDir inFile '.txt']);

% Psi-tet uses SI units so I will convert IDS data (in km/s and ms) to SI
% (m/s and s)
% dat.time = 1e-3 * dat.time;
% dat.vel = shiftVel + 1e3 * dat.vel;

% Channel Ranges Defined
switch configt
    case 1
        chan_ranget = [7:29]; % toroidal, 71 degree port
    case 2
        chan_ranget = [8:24]; % toroidal, mohawk port in midplane
    case 3
        chan_ranget = [10:27]; % toroidal, axial port at 135 degrees
    case 4
        chan_ranget = [8:28]; % toroidal, mohawk port perp.
end
switch configp
    case {1, 2, 3}
        chan_rangep = [46:63]; % 1: poloidal at toroidal angle 45 degrees
                               % 2: poloidal at toroidal angle 270 degrees
                               % 3: poloidal at toroidal angle 135 degrees
end

if tor % set channel range
    chan_range = chan_ranget;
    config = configt;
    tp = 't';
else
    chan_range = chan_rangep;
    config = configp;
    tp = 'p';
end

% Trim all Data for channel range
dat = trimRange(dat, chan_range);

% Turn all NaN into epsilon (setting everything to exactly zero throws an
% error in VisIt)
% eps = 1; % [m/s]
% nanVec = isnan(dat.vel(:));
% nanInd = find(nanVec);
% dat.vel(nanInd) = eps*ones(sum(nanVec), 1);

% Find coordinates and normalized vectors of IDS data

[pts, vec] = findIDSforVisIt(tor, config, dat);


%% Write Data to File

% for n = 1:length(dat.time)
    
    filename = [outDir inFile '.vtk'];
    fid = fopen(filename, 'w');
    
    fprintf(fid, '# vtk DataFile Version 3.0\n'); % (line 1)
    fprintf(fid, '%s\n', comment); % Comment (line 2)
    fprintf(fid, 'ASCII\n'); % (line 3)
    fprintf(fid, 'DATASET STRUCTURED_GRID\n'); % (line 4)
    fprintf(fid, ['DIMENSIONS 1 1 ' int2str(size(data, 1)) '\n']); % (line 5)
%     fprintf(fid, 'FIELD FieldData 1\n'); % this and next couple lines display time (line 5b)
%     fprintf(fid, 'TIME 1 1 double\n'); % (line 5c)
%     fprintf(fid, [num2str(dat.time(n), 4) '\n']); % (line 5d)
    fprintf(fid, ['POINTS ' int2str(size(data, 1)) ' FLOAT\n']); % (line 6)
    
    % Location Data
    for m = 1:size(pts, 1) % loop over all points
        fprintf(fid, [num2str(pts(m, :)) '\n']);
    end
    
    fprintf(fid, ' \n'); % (line 7)
    fprintf(fid, ['POINT_DATA ' int2str(size(data, 1)) '\n']); % (line 8)
    fprintf(fid, 'VECTORS V FLOAT\n'); % (line 9)
    
    % Vector Velocity Data
    for m = 1:size(vec, 1) % loop over all vectors
        scalarData = data(m, :) * vec(m, :)'; % dot 3D velocity vectors into chord = v dl
        fprintf(fid, [num2str(scalarData * vec(m, :)) '\n']);
    end
    
    fprintf(fid, ' '); % (line 10) last line in file
    
    fclose(fid); % close file
% end
    
    
    
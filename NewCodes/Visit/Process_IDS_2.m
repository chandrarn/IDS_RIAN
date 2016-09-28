%% Process_IDS_2
% Aaron Hossack
% Dec. 4th, 2013
%
% Update 5/28/14 - Due to positive velocity now defined as AWAY from the
% detector, minus sign added to velocity dotted with vec.
%
% Update 3/14/14 - The '2' version attempts to make a grid instead of
% discrete points.
% 
% First attempt at putting IDS data into .vtk format for VisIt.
% This version makes discrete vectors at discrete points
%
clear all; close all; clc;
addpath('~/IDS/Matlab/');
addAllThePaths;

%% Settings

shot = 12953010; % shot number corresponding to 'dat' file
shiftVel = 0; % [m/s] velocity shift
                  % 0 for 129793 and similar
                  % -7000 for 129499 and similar
                  % 16000 for 129591 and similar
tor = 0; % 1 = toroidal array (1:36), 0 = poloidal array (37:72)
outFile = '~/IDS/Visit/IDSvtk/ids';
    % shot number, tor/pol, and time points will be added to the file name
str = date;
comment = ['BD, ' str];

configt = 3; % 1 toroidal, 71 degree port
             % 2 toroidal, mohawk port in midplane
             % 3 toroidal, axial port at 135 degrees
             % 4 toroidal, mohawk port perp.
configp = 2; % 1 axial fiber @ toroidal angle 45
             % 2 axial fiber @ toroidal angle 270
             % 3 axial fiber @ toroidal angle 135

ptShift = 0.005; % [m], y-shift in arrows from reality so vectors don't overlay in VisIt
gridShift = 0.005; % [m], amount to shift cell points to either side of real points, ie: half cell width
lengthBeyond = 0.08; % [m], distance along chord beyond ring to display vector origin

%% Load and Prepare Data

load(['dat' num2str(shot)]);

% Psi-tet uses SI units so I will convert IDS data (in km/s and ms) to SI
% (m/s and s)
dat.time = 1e-3 * dat.time;
dat.vel = shiftVel + 1e3 * dat.vel;

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
eps = 1; % [m/s]
nanVec = isnan(dat.vel(:));
nanInd = find(nanVec);
dat.vel(nanInd) = eps * ones(sum(nanVec), 1);

%% Find coordinates and normalized vectors of IDS data

[pts, vec, ~] = findIDSforVisIt(tor, config, dat, ptShift, lengthBeyond);

%% Calculate Cells for VTK Unstructured Grid Format

[pts, vec] = pts2grid(pts, vec, tor, config, gridShift, lengthBeyond);

n_pts = size(pts, 1) / 2; % number of actual points
n_cells = n_pts - 1; % number of cells
n_vert = 4; % number of vertices in cells
sz_cells = n_vert + 1; % number of vertices + 1
cell_type = 9; % 'VTK_QUAD' = 9, from reference manual
v_eps = 1.0001; 

%% Write Data to File

for n = 1:length(dat.time)
    
    filename = [outFile '_' num2str(shot) tp '_' num2str(n) '.vtk'];
    fid = fopen(filename, 'w');
    
    fprintf(fid, '# vtk DataFile Version 2.0\n'); % (line 1)
    fprintf(fid, '%s\n', comment); % Comment (line 2)
    fprintf(fid, 'ASCII\n'); % (line 3)
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n'); % (line 4)
%     fprintf(fid, ['DIMENSIONS 1 1 ' int2str(size(dat.vel, 2)) '\n']); % (line 5)
    fprintf(fid, 'FIELD FieldData 1\n'); % this and next couple lines display time (line 5b)
    fprintf(fid, 'TIME 1 1 double\n'); % (line 5c)
    fprintf(fid, [num2str(dat.time(n), 4) '\n']); % (line 5d)
    fprintf(fid, ['POINTS ' int2str(size(pts, 1)) ' FLOAT\n']); % (line 6)
    
    % Location Data
    for m = 1:size(pts, 1) % loop over all points
        fprintf(fid, [num2str(pts(m, :)) '\n']);
    end
    fprintf(fid, ' \n');
    
    % Cells
    fprintf(fid, ['CELLS ' int2str(n_cells) ' ' int2str(sz_cells * n_cells) '\n']);
    for m = 1:n_cells
        fprintf(fid, [int2str([n_vert, m-1, m, m + n_pts, m + n_pts - 1]) '\n']); % index from zero!
    end
    fprintf(fid, ' \n');
    
    % Cell Types
    fprintf(fid, ['CELL_TYPES ' int2str(n_cells) '\n']);
    for m = 1:n_cells
        fprintf(fid, [int2str(cell_type) '\n']);
    end
    fprintf(fid, ' \n');
    
    % Vectors
%     fprintf(fid, ' \n'); % (line 7)
    fprintf(fid, ['POINT_DATA ' int2str(size(vec, 1)) '\n']); % (line 8)
    fprintf(fid, 'VECTORS V FLOAT\n'); % (li:ne 9)
    
    % Vector Velocity Data
    for m = 1:n_pts % loop over first half of vectors
        for p = 1:3
            fprintf(fid, [num2str(-dat.vel(n, m) * vec(m, p)) ' ']);
        end
        fprintf(fid, '\n');
    end
    for m = 1:n_pts % loop over second half of vectors
        for p = 1:3 % x, y, z components
            fprintf(fid, [num2str(-v_eps * dat.vel(n, m) * vec(m, p)) ' ']);
        end
        fprintf(fid, '\n');
    end
    
    fprintf(fid, ' '); % (line 10) last line in file
    
    fclose(fid); % close file
end
    
    
    
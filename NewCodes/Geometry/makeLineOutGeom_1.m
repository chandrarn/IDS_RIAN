%% makeLineOutGeom_1
% Aaron Hossack
% Dec. 4th, 2013
%
% Mofified to 'makeLineOutGeom_1.m' April 24th, 2014.
% This uses the same subfunction, 'findIDSforVisIt.m' as 'Process_IDS_2.m',
% but this script saves the lens location and end points to be used by
% VisIt and python to extract LineOut() data along simulated chords through
% simulations.
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

% shot = 12979310; % shot number corresponding to 'dat' file
% shiftVel = 0; % [m/s] velocity shift

tor = 1; % 1 = toroidal array (1:36), 0 = poloidal array (37:72)
outFile = '~/IDS/Geometry/coords';
    % tor(1)/pol(0), and config will be added to the file name

configt = 1; % 1 toroidal, 71 degree port
             % 2 toroidal, mohawk port in midplane
             % 3 toroidal, axial port at 135 degrees
             % 4 toroidal, mohawk port perp.
configp = 1; % 1 axial fiber @ toroidal angle 45
             % 2 axial fiber @ toroidal angle 270
             % 3 axial fiber @ toroidal angle 135

lengthBeyond = 0.01; % [m], distance along chord beyond ring to display vector origin
ptShift = 0; % not used for this version

%% Load and Prepare Data

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

dat.peaks = chan_range;

%% Find coordinates and normalized vectors of IDS chords

[pts, vec, origin] = findIDSforVisIt(tor, config, dat, ptShift, lengthBeyond);

%% Write Data to File

disp(['saving ' outFile int2str(tor) int2str(config) '.mat']);
save([outFile int2str(tor) int2str(config)], 'origin', 'pts', 'chan_range');

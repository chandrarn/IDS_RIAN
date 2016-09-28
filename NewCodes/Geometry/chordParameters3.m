% Aaron Hossack
% Jan 17th 2012
%
% Updated June 2013
clear all; close all; clc;
addpath('~/IDS/Matlab/');
%addAllThePaths;
cd('T:\RChandra\NewCodes\Geometry');
%
% This code takes calibration data from backlighting the new IDS fibers
% onto a screen to determine the angle, volume, and impact parameters of
% the various chords.
% clear all; close all; clc;
%% Settings
saveLength71 = 0; % save chrod lengths for 71 degree port chords

save1 = 0; % save 'impacts' as 'impacts1.mat'.  This configuration has the 
           % 'Upper/Long' fiber (smaller channel numbers) in the "old"
           % reentrant port at 71 degrees in machine coordinates.  Channel
           % 1 is oriented toward the top of the machine and channel 36 the
           % bottom.  The 'lower/short' fiber is on an axial port at 45
           % degrees on the X side of the machine.  Channel 37 is looking
           % toward the the small cones (inboard) and channel 72 is looking
           % toward the ring (outboard).
           
saveNim1 = 0; % save 'IDScoords.mat' with 'origin' and 'IDSchords' variables
           
save2 = 0; % save 'impacts' as 'impacts2.mat'.  This configuration has the 
           % 'Upper/Long' fiber in the upper mohawk port #6 with the fiber
           % fan in midplane.  Channel 1 is aimed down toward the bottom of
           % the machine and channel 36 is aimed up.  The 'lower/short'
           % fiber is on an axial port (either 45 or 270 degrees) with
           % channel 37 always aiming inboard, same as 'impacts1'.
           
saveLengthMohawk = 0; % save chord lengths for the above configuration
           
save3 = 0; % save 'impacts' as 'impacts3.mat'.  This is a special
           % configuration for velocity zeroing where the 'Upper/Long'
           % fiber is on the X side 135 degree port and the 'lower/short'
           % fiber is on the Y side 135 degree port.  Channels 1 and 37
           % both aim toward the geometric axis and 36 and 72 aim toward
           % the ring.
           
save4 = 0; % save 'impacts' as 'impacts4.mat'.  This has the long fiber in 
           % the mohawk port, just as in configuration 2, but rotated 90
           % degrees counter clockwise.  Channel 1 faces the y side
           % (positive z direction).

save5 = 1; % save 'impacts' as 'impacts5'. Upper mowhawk #6, lower mowhawk 
           % #6.
% raw = [247.5, 250.5, 253.5, 256.5, 259.5, 263.0,... % real unaltered data
%     265.5, 268.0, 271.0, 274.0, 277.0,...
%     280.0, 282.5, 285.0, 287.5, 290.5,...
%     293.0, 295.5, 298.5, 301.0, 304.0,...
%     306.5, 309.5, 312.0, 315.0, 318.0,...
%     321.0, 324.0, 327.0, 330.0, 332.5,...
%     336.0, 339.0];

% raw = [244.6, 247.5, 250.5, 253.5, 256.5, 259.5, 263.0,... % filling in extra 3 points
%     265.5, 268.0, 271.0, 274.0, 277.0,...
%     280.0, 282.5, 285.0, 287.5, 290.5,...
%     293.0, 295.5, 298.5, 301.0, 304.0,...
%     306.5, 309.5, 312.0, 315.0, 318.0,...
%     321.0, 324.0, 327.0, 330.0, 332.5,...
%     336.0, 339.0, 341.9, 344.8];
raw = 2.95 * [1:36];

% Center on zero degrees

angles = raw - raw(1) - 0.5*(raw(end) - raw(1)); % chord angles out of lens

R = 55.5; % cm, radius of machine

alpha = 180 - 2*angles; % angle from one end of chord, to geometric axis, to other end of chord
                        % (because the angles in a triangle must add up to
                        % 180)

r = 2 * R * sind(0.5 * alpha); % chord lengths

if saveLength71
    chordLength = r;
    save chordLength71.mat chordLength
end

x = r .* cosd(angles); % projection of chords on midplane diameter
% figure(14)
% plot(x)

% x = 222./(2.*(tand(angles).^2 + 1)); % projection of chords on midplane diameter

% r = x./cosd(angles); % length of chords

ratio = 0.02667/1.11; % radius to height ratio of cone

vol = pi*ratio^2*(1/3)*r.^3;

impacts1 = -R*sind(angles); % sign flip simply for orientation w.r.t machine
% figure(15)
% plot(impacts1)

%% For NIMROD comparison, find angles

y = r .* sind(angles); % dist. from diameter perp. up to chord end

theta = asind(y ./ R); % find angles w.r.t geom. axis from mid. of array

theta(1:3) = theta(1:3) + 2 .* (-90 - theta(1:3)); % mirror these angles past 90
theta(end-2:end) = theta(end-2:end) + 2 .* (90 - theta(end-2:end));

theta = theta + 251; % translate to machine coordinates

Rnim = 53.385; % NIMROD HIT-SI major radius

origin = [Rnim, 0, 71]; % coordinates of IDS wide angle lens

IDSchords = ones(length(theta), 3); % initialize array

IDSchords(:, 1) = Rnim * IDSchords(:, 1); % plug in R coordinates

IDSchords(:, 2) = 0 * IDSchords(:, 2); % zero out Z coordinate

IDSchords(:, 3) = theta' .* IDSchords(:, 3); % set toroidal angles

figure(12)
plot(IDSchords(:, 3), '+');
title('NIMROD coordinates 1, toroidal angles');

if saveNim1
    save IDScoords.mat origin IDSchords;
end

%% For axial port, find major radius where chord intersects midplane
% rPort = 36.8; % cm, approximate radius of axial port
% 
% rMAxis = 34; % cm, radius of magnetic axis
% 
% Rpol = 31.26; % cm, distance from midplane to wall
% 
% a0 = atand((rPort - rMAxis) / Rpol); % angle offset, ie: angle from z-hat to mag. axis
% 
% hypotn = hypot((rPort - rMAxis), Rpol); % distance directly from port to mag. axis
% 
% pAngles = angles + a0; % poloidal angles w.r.t. mag. axis
% 
% impacts1 = [impacts1, hypotn*sind(pAngles)];
rPort = 36.8; % cm, major radius location of port

zPort = 31.26 + 2.54 + 1.27; % cm, distance from magnetic axis to wall + thickness of window + thickness of wall

yPol = zPort * tand(angles(end:-1:1));

Rpol = rPort - yPol; % convert to major radius

impacts1 = [impacts1, Rpol];

figure(1)
plot(impacts1, '+');
title('Impacts 1');

if save1
    impacts = impacts1;
    save impacts1.mat impacts
end

%% Update - View with fan in midplane, mohawk port #6

theta_port = 31.87; % degrees, port axis w.r.t machine radius

theta_chords = angles + theta_port;

y_tor = R * sind(theta_chords);

impacts2 = impacts1; % copy over previous impacts array

impacts2(1:36) = y_tor;

figure(2)
plot(impacts2, '+');
title('Impacts 2');

if save2
    impacts = impacts2;
    save impacts2.mat impacts
end

%% Update 3-18-15 - Chord lengths for the above configuration

alpha2 = 180 - 2*theta_chords; % angle from one end of chord, to geometric axis, to other end of chord
                        % (because the angles in a triangle must add up to
                        % 180)

r2 = 2 * R * sind(0.5 * alpha2); % chord lengths                        
                        
if saveLengthMohawk
    chordLength = r2;
    save chordLengthMohawk.mat chordLength;
end

%% Velocity Zeroing configuration

impacts3 = [impacts1(37:72), impacts1(37:72)];

figure(3)
plot(impacts3, '+');
title('Impacts 3');

if save3
    impacts = impacts3;
    save impacts3.mat impacts
end

%% Mohawk #6, chords perpendicular to midplane

half_ch_length = R * cosd(theta_port); % half chord length of axis of lens

displacement = half_ch_length * sind(angles); % displacement axially from magnetic axis.

impacts4 = [displacement(end:-1:1), impacts2(37:72)];


if save4
    figure(4)
    plot(impacts4, '+');
    title('Impacts 4');

    impacts = impacts4;
    save impacts4.mat impacts
end

%% Top Mohawk #6, Bottom Mohawk #6

theta_port = 31.87; % degrees, port axis w.r.t machine radius

theta_chords = -angles + theta_port;

y_tor = R * sind(theta_chords);

impacts5(1:36) = y_tor;
theta_port = -31.87; % degrees, port axis w.r.t machine radius

theta_chords = angles + theta_port;

y_tor = R * sind(theta_chords);

impacts5(37:72) = y_tor;

if save5
    figure(4)
    plot(impacts5, '+');
    title('Impacts 5');

    impacts = impacts5;
    save impacts5.mat impacts
end





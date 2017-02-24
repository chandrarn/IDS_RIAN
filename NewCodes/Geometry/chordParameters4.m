% Aaron Hossack
% Jan 17th 2012
%
% Updated June 2013
%clear all; 
close all; clc;
addpath('~/IDS/Matlab/');
%addAllThePaths;
%cd('T:\RChandra\NewCodes\Geometry');
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
           
saveNim1 = 1; % save 'IDScoords.mat' with 'origin' and 'IDSchords' variables
           
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

save5 = 0; % save 'impacts' as 'impacts5'. Upper mowhawk #6, lower mowhawk 
           % #6.

save6 = 0; % save 'impacts' as 'impacts6', Upper mohawk #1-#10, with secant
           % and nimrod coordinates

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

% assuming IDS is in port 6 and -6
CoordRot = 225;
theta_port = 31.87; % degrees, port axis w.r.t machine radius
anglesTop = angles + theta_port;
anglesBot = angles - theta_port;
terminalAng1 = 180-2.*anglesTop + 59;
terminalAng2 = 180-2.*anglesBot - 59;
theta1 = terminalAng1;
theta2 = terminalAng2;
theta = [terminalAng1(end:-1:1), terminalAng2];
  theta = mod(theta +CoordRot ,360); % translate to machine coordinates
% 

% assuming IDS is in a radial port
% y = r .* sind(angles ); % dist. from diameter perp. up to chord end
% 
% thetaR = asind(y ./ R); % find angles w.r.t geom. axis from mid. of array
% 
% thetaR(1:3) = thetaR(1:3) + 2 .* (-90 - thetaR(1:3)); % mirror these angles past 90
% thetaR(end-2:end) = thetaR(end-2:end) + 2 .* (90 - thetaR(end-2:end));
% originR = [53.385, 1, 71]

% 
thetaR = 180-2*angles +71;
 originR = [53.385, 1, 71]

 Rnim = 53.385; % NIMROD HIT-SI major radius
% 
 origin = [Rnim, 0, 71+225 ; Rnim, 0, -59]; % coordinates of IDS wide angle lens
origin(1:36,1) = ones(1,36).*Rnim;
origin(1:36,2) = ones(1,36).*1;
origin(1:36,3) = ones(1,36).*(59+CoordRot);
origin(37:72,1) = ones(1,36).*Rnim;
origin(37:72,2) = ones(1,36).*1;
origin(37:72,3) = ones(1,36).*(-59+CoordRot);

origin(:,1:2) = origin(:,1:2).*.01;

IDSchords = ones(length(theta), 3); % initialize array

IDSchords(:, 1) = Rnim * IDSchords(:, 1) .* .01; % plug in R coordinates

IDSchords(:, 2) = 1 * IDSchords(:, 2) .* .01; % zero out Z coordinate

IDSchords(:, 3) = theta' .* IDSchords(:, 3); % set toroidal angles




figure(12)
polar([ones(1,length(theta1)).*originR(3); thetaR].*(pi/180),[ones(36,1).*Rnim ones(36,1).*Rnim]' ,'-*');
[ones(1,length(theta1)).*originR(3); thetaR].*(pi/180),[ones(36,1).*Rnim ones(36,1).*Rnim]
title('71^o Port')
% plot(IDSchords(:, 3), '+');
figure(13)
polar([ones(1,length(theta1)).*origin(1,3); theta(1:36)].*(pi/180),[ones(36,1).*Rnim ones(36,1).*Rnim]' ,'-*');
[ones(1,length(theta1)).*origin(1,3); theta(1:36)].*(pi/180),[ones(36,1).*Rnim ones(36,1).*Rnim]
title('IDS upper/lower fibers, Lab Frame');
hold on;
polar([ones(1,length(theta2)).*origin(37,3); theta(37:end)].*(pi/180),[ones(36,1).*Rnim ones(36,1).*Rnim]' ,'-*');
[ones(1,length(theta2)).*origin(37,3); theta(37:end)].*(pi/180),[ones(36,1).*Rnim ones(36,1).*Rnim]

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


if save1
    figure(1)
plot(impacts1, '+');
title('Impacts 1');

    impacts = impacts1;
    save impacts1.mat impacts
end

%% Update - View with fan in midplane, mohawk port #6

theta_port = 31.87; % degrees, port axis w.r.t machine radius

theta_chords = angles + theta_port;

y_tor = R * sind(theta_chords);

impacts2 = impacts1; % copy over previous impacts array

impacts2(1:36) = y_tor;



if save2
    figure(2)
plot(impacts2, '+');
title('Impacts 2');
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



if save3
    figure(3)
plot(impacts3, '+');
title('Impacts 3');
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
    figure(5)
    plot(impacts5, '+');
    title('Impacts 5');

    impacts = impacts5;
    save impacts5.mat impacts
end

%% all the bottom mohawk
if save6 

Origin(1:10,1) = ((13.736) + (0:9)'.*9) ;% normalize later + 225; % ports, in Normal Machine Coordinates (Tm)
Gamma(1:10,1) = 180 - Origin(1:10,1); % angle between origin-port and origin-convergence point
a = 57.785; %Radial location that the ports converge to 
%R = 60.96; % R machine ( MAY BE TANK R, NOT PLASMASPRAY R
R = 55.5; % 
c(1:10,1) = sqrt( R.^2 + a.^2 -2*R*a*cosd( Gamma )); % chord from port to convergence point 
B(1:10,1) = asind(a*sind(Gamma)./c); % angle off of the radial that the port is
% angle between origin-port and port-convergence point
GammaB = 180-2.*B; % angle between origin-port and origin-(intersection of port-convergence point
% with wall of chamber)
Lengths(1:10,1) = 2* R*sind(GammaB); % apply law of sines for a triangle subtending a circle,
% get length of secant across circle from port to wall
Origin(1:10,2) = Origin(1:10)+GammaB; % ending theta locations, still John coordinates
figure(7);
%Machine coords

Origin = mod(Origin +CoordRot,360);
OutputI = [ ones(3,1).* Rnim .* .01 ones(3,1) .* .01 Origin(1:3,1)];
OutputF = [ ones(3,1).* Rnim .* .01 ones(3,1) .* .01 Origin(1:3,2)];

IDSchords = [IDSchords ; OutputF];
origin = [origin ; OutputI];

h=polar(Origin'.*(pi/180),[ones(10,1).*55.785 ones(10,1).*57]' ,'-*');

title('Mohawk lines of sight, Lab frame');

save coords06.mat IDSchords origin
end




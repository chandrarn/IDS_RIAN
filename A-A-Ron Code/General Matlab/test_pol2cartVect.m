% Aaron Hossack
% 
% Test vector pol2cart converter.

clear all; close all; clc;

theta = [pi/4 3/4*pi 5/4*pi -pi/4];
rho = [2 2 2 2];

Vt = [-1 -1 -1 -1];
Vr = [-1 -1 -1 -1];

[X, Y, Vx, Vy] = pol2cartVect(theta, rho, Vt, Vr);

quiver(X, Y, Vx, Vy);
hold on
grid on
set(gca, 'XLim', [-5 5], 'YLim', [-5 5]);
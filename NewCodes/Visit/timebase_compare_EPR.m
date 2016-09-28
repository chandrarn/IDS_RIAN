% Aaron Hossack
clear all; close all; clc;
addpath('~/IDS/Matlab/');
addAllThePaths;

nimTimeFile = '/home/aaron/IDS/Visit/NIMvtk';
nimInjFile = '/home/aaron/IDS/NimrodData/fin_beta_ramp/24x24_14_injRamp.mat';

tetTimeFile = '/home/aaron/IDS/Visit/TETvtk';

% Generate PSI-TET injector current

f=14.5; % [kHz]
amp = 21; % [kA]
tet.injt = linspace(0, 0.6, 10000);
tet.injx = amp * sin(2*pi*f*tet.injt);


%% Load and Plot PSI-TET data times
load([tetTimeFile '/t7129499']);
tet.t = t;
tet.x = ones(1, length(t));

figure(1)
plot(tet.t, tet.x, '+b');
hold on;
plot(tet.injt, tet.injx, 'b');
title('PSI-TET');
grid on;
xlabel('time [ms]');

%% NIMROD

load([nimTimeFile '/t8129499']);
nim.t = t;
nim.x = ones(1, length(t));

load(nimInjFile);
nim.injx = -nimsave.inj.I1;
nim.injt = 1e3 * nimsave.time;
nim.injy = nimsave.inj.I2;

figure(2)
plot(nim.t, nim.x, 'or');
hold on;
plot(nim.injt, nim.injx, '*-r');
plot(nim.injt, nim.injy, '*-b');
title('NIMROD');
legend('dump file times', 'X inj', 'Y inj');
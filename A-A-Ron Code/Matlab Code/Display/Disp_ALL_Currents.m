% Aaron Hossack
% December 31st, 2013
%
% Script for displaying NIMROD currents to determine phasing
clear all; close all; clc;

%% IDS Settings
shot = '12949910';

cd('S:\MC_IDS\Matlab Code\Data Repository');
load(['dat' shot]);
IDStime = dat.time;
IDSItorTime = dat.ItorTime;
IDSiinjxTime = dat.iinjxTime;
IDSiinjx = dat.iinjx;
IDSItor = dat.Itor;

%% NIMROD Settings
matFile = 'NIM_IDS'; % file name
fileID = '13080800'; % folder where mat file is stored

cd(['S:\MC_IDS\Matlab Code\NIMROD\' fileID]);
load(matFile);
NIMtime = 1e3 * time;
NIMiinjx = 1e-3 * iinjx;
NIMItor = 1e-3 * Itor;

%% PSI-TET Settings
psiFile = 'PSITET_eta16_Hall';
dataFolder = 'aaronData_131029';

cd(['S:\MC_IDS\Psi TET Comparisons\' dataFolder]);
load(psiFile); % contains toroidal current
TETtime = 1e3 * Time;
TETItor = 1e-3 * Tcurr; % toroidal current

% Create Injector Current
f = 14.5; % [Hz]
w = 2*pi*f; % omega
amp = 20; % [kA] amplitude for plotting
TETiinjx = amp * sin(w * TETtime);

% Extract time points of Psi-Tet data
% n_time = 1; % number of time points
n_time = 40; % number of time points
TETdataTime = zeros(1, length(n_time));
for m = 1:n_time
        if m < 10
            setname = ['out_000' int2str(m) '.xmf'];
        elseif m < 100
            setname = ['out_00' int2str(m) '.xmf'];
        end
        
        text = fileread(setname);
        tmp = regexp(regexp(text,'Time Value="[0-9.]*"','match'),'[0-9.]*','match');
        TETdataTime(m) = str2num(char(tmp{1}));
        clear text tmp;
end
TETdataTime = 1e3 * TETdataTime; % convert from [s] to [ms]

%% Plotting
S = get(0, 'ScreenSize');
h1 = figure('Visible', 'on', 'Name', 'All Currents', 'Position',...
    [0.05*S(3), 0.1*S(4), 0.9*S(3), 0.8*S(4)], 'Color', [1 1 1]);

%% IDS
ax1(1) = subplot(2,3,1);
plot(IDSItorTime, IDSItor, '-*b');
hold on;
grid on;
plot(IDSiinjxTime, IDSiinjx, '-*b');
title('IDS');
xlabel('time [ms]');
ylabel('Current [kA]');

ax1(2) = subplot(2,3,4);
plot(IDStime, 1:length(IDStime), '-*b')
hold on;
grid on;
xlabel('time [ms]');
ylabel('Index');

linkaxes(ax1, 'x');

%% NIMROD
ax2(1) = subplot(2,3,2);
plot(NIMtime, NIMItor, '-*b');
hold on;
grid on;
plot(NIMtime, NIMiinjx, '-*b');
title('NIMROD');
xlabel('time [ms]');
ylabel('Current [kA]');

ax2(2) = subplot(2,3,5);
plot(NIMtime, 1:length(NIMtime), '-*b')
hold on;
grid on;
xlabel('time [ms]');
ylabel('Index');

linkaxes(ax2, 'x');

%% PSI-TET
ax3(1) = subplot(2,3,3);
plot(TETtime, TETItor, '-*b');
hold on;
grid on;
plot(TETtime, TETiinjx, '-*b');
title('PSI-TET');
xlabel('time [ms]');
ylabel('Current [kA]');

ax3(2) = subplot(2,3,6);
plot(TETtime, 1:length(TETtime), '-*b')
hold on;
grid on;
plot(TETdataTime, 1:length(TETdataTime), '-*r')
xlabel('time [ms]');
ylabel('Index');

linkaxes(ax3, 'x');
% Aaron Hossack
% December 31st, 2013
%
% Script for displaying NIMROD currents to determind phasing
clear all; close all; clc;

matFile = 'NIM_IDS'; % file name
fileID = '13080800'; % folder where mat file is stored



cd(['S:\MC_IDS\Matlab Code\NIMROD\' fileID]);
load(matFile);

figure(1)

ax(1) = subplot(2,1,1);
plot(Itor, '-*b');
hold on;
plot(iinjx, '-*b');

ax(2) = subplot(2,1,2);
plot(time, '-*b')
hold on;

linkaxes(ax, 'x');
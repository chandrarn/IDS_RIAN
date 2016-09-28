% Aaron Hossack

clear all; close all; clc;
addpath('~/IDS/Matlab/');
addAllThePaths;

sim = 'nim'; % 'tet' or 'nim'
tor = 1; % tor = 1, pol = 0
config = 2; % config of the primary dataset


%% LOAD LINEOUTS

if strcmp(sim, 'tet')
    simFolder = 'TETvtk';
elseif strcmp(sim, 'nim')
    simFolder = 'NIMvtk';
end
cd(['/home/aaron/IDS/Visit/' simFolder]);

matList = dir(['LO' sim int2str(tor) '_' int2str(config) '_*.mat']); % all times for this fiber and configuration

nTimes = size(matList, 1);

for m = 1:nTimes
    load(matList(m).name); % loads 'lineVals' and 'chan_range'
    [~, ~, timeOld] = strread(matList(m).name, '%s %s %d', 'delimiter', '_');
    if timeOld < 10
        timeNew = ['000' num2str(timeOld)];
    elseif timeOld < 100
        timeNew = ['00' num2str(timeOld)];
    elseif timeOld < 1000
        timeNew = ['0' num2str(timeOld)];
    else
        timeNew = num2str(timeOld);
    end
    cd('revised');
    save(['LO' sim int2str(tor) '_' int2str(config) '_' timeNew '.mat'], 'lineVals', 'chan_range');
    cd('..');
    
end

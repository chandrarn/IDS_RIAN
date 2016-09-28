% Aaron Hossack
% Sept. 26th, 2014
% 
% Correct dat files which have the wrong 'impacts' file
clear all; close all; clc;
cd('T:\IDS\Data Repository');

%% Settings

refShot = -1; % shot with correct impacts file in tree

% possible shots that have already been corrected. It's ok if they have not
datShots = [140923006:140923026, 140924005:140924038, 140925005:140925034];
% datShots = 140924033;

%% Load Impacts
import MDSplus.*;

HitTree = Tree('hitsi3', refShot);
Data = HitTree.getNode('\IDS_IMPACTS');
impacts = NATIVEvalue(Data.getData());
mdsclose();

for n = 1:length(datShots)
    try
        load(['dat' num2str(datShots(n))]);
        
        dat.impacts = impacts;
        save(['dat' num2str(datShots(n))], 'dat');
        disp(['saved ' num2str(datShots(n))]);
        clear dat;
    end
end
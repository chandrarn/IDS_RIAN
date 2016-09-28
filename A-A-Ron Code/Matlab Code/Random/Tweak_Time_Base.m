% Correct Time Base
%
% This is a single-use script to slightly adjust the time bases of all
% previously corrected shots.  The Phantom camera places the time zero at
% the END of the exposure, so this moves the time base back by 1/2 the
% exposure so the data point is centered.

clear all; close all; clc;
cd('S:\MC_IDS\Matlab Code\Data Repository');

shotStart = 13060501;
shotEnd = 13071902;

special = [129499, 130710001, 130710004];

for n = [shotStart:shotEnd, special]
    try
        load(['dat' num2str(n)]);
        n
        cd('T:\PhantomMovies');
        dat.shotRef
        load(['tag' num2str(dat.shotRef)]);
        d=1
        cd('S:\MC_IDS\Matlab Code\Data Repository');
        e=1
        
        dat.time = dat.time - 0.5*exposure;
        f=1
        0.5*exposure
        
%         save(['dat' num2str(n)], 'dat');
        clear corrected exposure dat;
        g=1
        disp(['Altered ' num2str(n)]);
        h=1
    end
end
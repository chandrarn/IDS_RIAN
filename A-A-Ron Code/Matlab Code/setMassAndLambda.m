% Simple script to set the ion mass and wavelengths in the tree
clear all; close all; clc;

saveToShot = [128457:128540, -1]; % shot to save settings to.  Use the latest version of 'BackUpdate' to save to more

ION_MASS = [12];
LINE_LAMBDA = 1e-9*[465.025, 465.025];
IDS_VOLTAGE = 900;

mdsconnect('landau.hit');
for n = 1:length(saveToShot)
    mdsopen('analysis', saveToShot(n));
    disp(['Saving to shot ' num2str(saveToShot(n))]);
    
%     mdsput('\ION_MASS', '$', ION_MASS);
    mdsput('\LINE_LAMBDA', '$', LINE_LAMBDA);
%     mdsput('\IDS_VOLTAGE', '$', IDS_VOLTAGE);
    mdsclose();
end
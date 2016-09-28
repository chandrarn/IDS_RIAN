% When the fast camera is removed from service on IDS, run this script to
% blow away nodes in -1 tree that are no longer relevent.
clear all; close all; clc;

shots = -1; % shots to clean out, usually just -1

for n = 1:length(shots)
    output = [];
    disp(['Editing ANALYSIS Tree for shot ' num2str(shots(n))]);
    
    % Open Tree for Edit
    output = [output, mdsvalue(['tcl("edit ANALYSIS/shot=' num2str(shots(n)) '")'])];
    
    % Navigate to SPECTROSCOPY node
    output = [output, mdsvalue('tcl("set def .spectroscopy.ids")')];
    
    % Nuke nodes 
    output = [output, mdsvalue('tcl("delete node IDS_PIX_SP /noconfirm")')];
    output = [output, mdsvalue('tcl("delete node IDS_FWHM /noconfirm")')];
    output = [output, mdsvalue('tcl("delete node IDS_CENTER /noconfirm")')];
    output = [output, mdsvalue('tcl("delete node IDS_BOUNDS /noconfirm")')];
    output = [output, mdsvalue('tcl("delete node IDS_PEAKS /noconfirm")')];
    
    % Write and close
    output = [output, mdsvalue('tcl("write")')];
    output = [output, mdsvalue('tcl("close")')];
    
end
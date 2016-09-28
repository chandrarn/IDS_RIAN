%% addAllThePaths
% Aaron Hossack
% March 13th, 2014
%
% This script adds many paths that may be required on the Linux system and
% any other system, like the 'T' drive.
%
warning('off', 'all'); % disable all warnings

addpath('/home/aaron/IDS/Matlab/'); % this folder
addpath('/home/aaron/IDS/Matlab/LM Method'); 
% addpath('~/IDS/Matlab/lsqcurvefit/'); % lsqcurvefit
addpath('/home/aaron/IDS/Matlab/extrema/'); % extrema
addpath('/home/aaron/IDS/Visit/'); % scripts related to Visit
addpath('/home/aaron/IDS/Fitting/'); % IDS basic data analysis/fitting
addpath('/home/aaron/IDS/Calibration/'); % calibration codes
addpath('/home/aaron/IDS/IDSdata/'); % where IDS data is stored
addpath('/home/aaron/IDS/Display/'); % display IDS data
addpath('/home/aaron/IDS/Geometry/'); % chord geometries
addpath('/media/alfventemp/PhantomMovies'); % all Phantom camera movies, Linux
addpath('T:\PhantomMovies'); % all Phantom camera movies, Windows
addpath('T:\RChandra\A-A-Ron Code\Matlab Code\Core Fitting Codes');
addpath('T:\RChandra\A-A-Ron Code\General Matlab');
addpath('T:\RChandra\A-A-Ron Code\General Matlab\LM Method');
addpath('T:\IDS\Data Analysis');
addpath('C:\Program Files\MDSplus\MATLAB');
addpath('T:\RChandra\A-A-Ron Code\Matlab Code\Core Fitting Codes\lsqcurvefit');

warning('on', 'all'); % enable all warnings
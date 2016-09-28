function AddAllThePaths
    if isunix % if we're on linux
        addpath('/media/alfventemp/IDS/General Matlab/extrema');
        addpath('/media/alfventemp/IDS/General Matlab/LM Method');
        addpath('/media/alfventemp/IDS/Calibration');
        addpath('/media/alfventemp/RChandra/A-A-Ron Code/Matlab Code/Core Fitting Codes/lsqcurvefit');
        addpath('/media/alfventemp/RChandra/A-A-Ron Code/Matlab Code/Core Fitting Codes');
        addpath('/media/alfventemp/PhantomMovies/');
    else
        addpath('T:\IDS\General Matlab\extrema');
        addpath('T:\IDS\General Matlab\LM Method');
        addpath('T:\IDS\Calibration');
        addpath('T:\RChandra\A-A-Ron Code\Matlab Code\Core Fitting Codes\lsqcurvefit');
        addpath('T:\RChandra\A-A-Ron Code\Matlab Code\Core Fitting Codes');
        addpath('T:\PhantomMovies\');
    end
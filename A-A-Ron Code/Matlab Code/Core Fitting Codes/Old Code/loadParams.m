function[param, options] = loadParams(shot, line)

param.shotRef = shot;

% Load Parameters from Tree
HitTree=Tree('hitsi',shot);

Node=HitTree.getNode('\IDS_PEAKS'); % Technically, this all can be done in one line.
para.peaks=NATIVEvalue(Node.getData());  % [pixels] (n_channels x 5) chan number, x, y, sig_x, sig_y
if length(line) == 2 % fit to a doublet
    param.dbl = 1; % flag to tell Gaussian fitting routine to fit doublet
else
    param.dbl = 0; % flag indicating single line- not a doublet
end

IonMass = (1.66e-27)*NATIVEvalue(HitTree.getNode('\IDS_MASS').getData()); % [AMU] => [kg]
param.IonMass = IonMass(line(1)); % select mass of ion of interest- for doublet, first one

param.CalLam = NATIVEvalue(HitTree.getNode('\CAL_LAMBDA').getData()); % [m]
LineLam = NATIVEvalue(HitTree.getNode('\IDS_LAMBDA').getData()); % [m]
param.LineLam = LineLam(line); % skip calibration lambda, which is first
% for doublet, this automatically selects both lines

param.PIX_SP = NATIVEvalue(HitTree.getNode('\IDS_PIX_SP').getData()); % [m] (1 x n_channels) wavelength spacing per channel

param.REL_INT = NATIVEvalue(HitTree.getNode('\IDS_REL_INT').getData()); % INVERSE relative intensity, 
                                          % ie: multiply by fit area to calibrate

% Adjust Center from calibration line to line of interest- NB: for DOUBLET,
% center is arbitrarily chosen to be the smaller wavelength line ***
param.Center = param.peaks(:, 3) - (param.CalLam - LineLam(line(end))) ./ param.PIX_SP;

% Establish Constants
param.limits = [0 100; -40e3 40e3];
    % limits for sanitizing final output- throws out data  that exceeds the
    % above limits for [temp; vel]
param.n_chan = size(param.peaks, 1);

param.calcError = 0; % calculate uncertainty by LM method - presently broken
param.ampThresh = 8; % [pixels], amplitude threshold
param.xTol = 0.5; % [pixels], tolerance of recentering on line in spatial dimension
param.yTol = 3.5; % [pixels], tolerance for velocity offset
param.xWing = 2; % real-space domain for fitting. 1 = 3 wide, 2 = 5 wide
param.yWing = 7; % wavelength-space domain. 5 = 11 wide, etc...
options = optimsetv61('lsqcurvefit'); % set options for curve fitting
param.kBoltz = 1.3807e-23; % Boltzmann's Constant [SI]
param.c = 2.99792458e8; % speed of light [SI]

end
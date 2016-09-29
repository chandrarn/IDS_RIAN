function[param, options] = loadParams(shot, line,hitsi3)

    param.shotRef = shot;

    import MDSplus.*;
    if(~hitsi3) % I ASSUME THAT ALL THE NODES I WANT EXIST ON HITSI AND HITSI3
        HitTree = Tree('hitsi',shot);
    else
        HitTree= Tree('hitsi3',shot);
    end
    
    
   Data = HitTree.getNode('\IDS_PEAKS');
    param.peaks = NATIVEvalue(Data.getData()); % [pixels] (n_channels x 5) chan number, x, y, sig_x, sig_y
    %param.peaks(9:27,2)=[41;44;47.5;50.5;54;57;61;64;67.5;71;74;78;81;88;91;95;99;102.5;106];
     
    if length(line) == 2 % fit to a doublet
        param.dbl = 1; % flag to tell Gaussian fitting routine to fit doublet
    else
        param.dbl = 0; % flag indicating single line- not a doublet
    end
    
    Data = HitTree.getNode('\IDS_MASS');
    IonMass = (1.66e-27)*NATIVEvalue(Data.getData()); % [AMU] => [kg]
    param.IonMass = IonMass(line(1)); % select mass of ion of interest- for doublet, first one
    
    Data = HitTree.getNode('\CAL_LAMBDA');
    param.CalLam = NATIVEvalue(Data.getData()); % [m]
    Data = HitTree.getNode('\IDS_LAMBDA');
    LineLam = NATIVEvalue(Data.getData()); % [m]
    param.LineLam = LineLam(line); % skip calibration lambda, which is first
    % for doublet, this automatically selects both lines

    Data = HitTree.getNode('\IDS_PIX_SP');
    param.PIX_SP = NATIVEvalue(Data.getData()); % [m] (1 x n_channels) wavelength spacing per channel
    
    Data = HitTree.getNode('\IDS_REL_INT');
    param.REL_INT = NATIVEvalue(Data.getData()); % INVERSE relative intensity, 
                                              % ie: multiply by fit area to calibrate

    % Adjust Center from calibration line to line of interest- NB: for DOUBLET,
    % center is arbitrarily chosen to be the smaller wavelength line ***
    param.Center = param.peaks(:, 3)+10 - (param.CalLam - LineLam(line(end))) ./ param.PIX_SP;

    % Establish Constants
    param.limits = [0 100; -40e3 40e3];
        % limits for sanitizing final output- throws out data  that exceeds the
        % above limits for [temp; vel]
    param.n_chan = size(param.peaks, 1);

    param.calcError = 1; % calculate uncertainty by LM method - presently broken
    param.ampThresh = 8; % [pixels], amplitude threshold
    param.xTol = 0.5; % [pixels], tolerance of recentering on line in spatial dimension
    param.yTol = 3.5; % [pixels], tolerance for velocity offset
    param.xWing = 2; % real-space domain for fitting. 1 = 3 wide, 2 = 5 wide
    param.yWing = 7; % wavelength-space domain. 5 = 11 wide, etc...
    options = optimsetv61('lsqcurvefit'); % set options for curve fitting
    param.kBoltz = 1.3807e-23; % Boltzmann's Constant [SI]
    param.c = 2.99792458e8; % speed of light [SI]
    
    
    mdsclose();

end
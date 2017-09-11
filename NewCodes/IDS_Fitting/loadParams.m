function[param, options] = loadParams(shot, line, hitsi3, useTree)
    % This code loads in calibration parameters from the tree into the
    % "param.XXXXX" structure. It further sets fitting parameters like the
    % minimum gaussian amplitude to be fit, and the window for BD filtering
    % and Gaussian fittings. It finally calculates the positions of the
    % lines to be analyzed, based on the position of line 3, the second OII
    % line.
    % It can load calibration data from outside the tree if necessary.
    
    
    %y0 offset, set the cal line the second O-II line (in between
    %C-Doublet)
    y0 = 1;11;8;
    x0 = 0;-4; % -2.5;
    if useTree % try using the tree, if that doesnt work, try importing callibration from file
        import MDSplus.*;
        Conn=Connection('landau.hit');
        if(~hitsi3) % I ASSUME THAT ALL THE NODES I WANT EXIST ON HITSI AND HITSI3
            Conn.openTree('hitsi', shot);
        else
            Conn.openTree('hitsi3', shot)
        end

        Data = Conn.get('\IDS_PEAKS');
        param.peaks = double(NATIVEvalue(Data.getFloatArray)); % [pixels] (n_channels x 5) chan number, x, y, sig_x, sig_y
        param.peaks = reshape(param.peaks,[],5); % sometimes this array is squeezed into a vector
        param.peaks(:,3) = param.peaks(:,3) + y0;
        param.peaks(:,2) = param.peaks(:,2) + x0;
        % Adjustments to PEAKS go here, for shifting X0 or Y0 positions/etc
        
        assignin('base','OrigPeaks',param.peaks);
        
        Data = Conn.get('\IDS_MASS');
        IonMass = (1.66e-27) * double(NATIVEvalue(Data.getFloatArray)); % [AMU] => [kg]
        param.IonMass = IonMass(line); % select mass of ion of interest- for doublet, first one

        Data = Conn.get('\CAL_LAMBDA');
        param.CalLam = double(NATIVEvalue(Data.getFloatArray)); % [m]
        Data = Conn.get('\IDS_LAMBDA');
        LineLam = double(NATIVEvalue(Data.getFloatArray)); % [m]
        param.LineLam = LineLam;
        % param.LineLam = LineLam(line); % This has to be put at the end,
        % we need all the lines to calculate position.
        
        Data = Conn.get('\IDS_PIX_SP');
        param.PIX_SP = double(NATIVEvalue(Data.getFloatArray)); % [m] (1 x n_channels) wavelength spacing per channel
       % Temporary adjustments to PIX_SP go here occationally
        
        Data = Conn.get('\IDS_REL_INT');
        param.REL_INT = double(NATIVEvalue(Data.getFloatArray)); % INVERSE relative intensity, 
        
        
                                       % ie: multiply by fit area to calibrate
        Data = Conn.get('\IDS_IMPACTS');
        param.impacts = double(NATIVEvalue(Data.getFloatArray));
        mdsclose();
    else
        %cd('T:\IDS\Calibration\Calibration repository')
%         temp = importdata(['T:\IDS\Calibration\Calibration repository\' int2str(shot) '.mat'] ); % load in params from file
        %temp = importdata(['C:\Users\Rian\Documents\MATLAB\thosematfilestho\Shot ' int2str(shot) '.mat'] ); % load in params from file
        temp = importdata(['dat' num2str(shot) '10.mat']);
        %assignin('base','temp',temp);
        param.peaks = temp(1).param.peaks;
        param.IonMass = temp(1).param.IonMass(line);
        param.CalLam = temp(1).param.CalLam;
        param.LineLam = [4.64181027837185e-07;4.64741788164247e-07;4.64913483710916e-07;4.65024612594789e-07];%temp(1).param.LineLam(line);
        LineLam=param.LineLam;
        param.PIX_SP = temp(1).param.PIX_SP;
        param.REL_INT = temp(1).param.REL_INT;
        param.Impacts = temp(1).param.impacts;
        clear temp
        assignin('base','param',param)
        
    end

%     if length(line) == 2 % fit to a doublet
%         param.dbl = 1; % flag to tell Gaussian fitting routine to fit doublet
%     else
        param.dbl = 0; % flag indicating single line- not a doublet
%     end

    %offset the line
    %param.peaks(:,3) = param.peaks(:,3) +y0;
    % slightly more clevel offset, takes into account curvature. Kinda.
    %param.CalLam = 464.913e-9 -(y0 - mean(param.peaks(:,3)))*mean(param.PIX_SP)
    
    

    param.shotRef = shot;
    % Adjust Center from calibration line to line of interest- NB: for DOUBLET,
    % center is arbitrarily chosen to be the smaller wavelength line ***
    % Ajust to center OII line
    for n = 1:length(line)
        %param.Center(:, n) = param.peaks(:, 3) + lineShift(n);
        for i = 1:size(param.peaks,1)
            if floor(log10(param.PIX_SP(i)))~=-11
                 temp = 1.18e-11;
            else
                temp = param.PIX_SP(i);
            end
            param.Center(i, n) = param.peaks(i, 3)+(param.LineLam(2) - param.LineLam(line(n))) ./ temp; %mean(param.PIX_SP);%LineLam replaces CalLam, Mean(pixsp) underpredicts
         end
    end
    
    param.LineLam = LineLam(line); % This needs to be after center finding
        
    % Establish Constants
    param.limits = [0 100; -40e3 40e3];
        % limits for sanitizing final output- throws out data  that exceeds the
        % above limits for [temp; vel]
    param.n_chan = size(param.peaks, 1);

    param.calcError = 1; % calculate uncertainty by LM method - presently broken
    param.ampThresh = 1; % [pixels], amplitude threshold TEMP SET TO ONE
    param.xTol = 0.5; % [pixels], tolerance of recentering on line in spatial dimension
    param.yTol = 3.5; % [pixels], tolerance for velocity offset
    param.xWing = 2; % real-space domain for fitting. 1 = 3 wide, 2 = 5 wide
    param.yWing = 7; % wavelength-space domain. 5 = 11 wide, etc...
    try
        options = optimsetv61('lsqcurvefit'); % set options for curve fitting
    catch
        options = optimset('lsqcurvefit'); % set options for curve fitting
    end
    % MODIFY MAXIMUM LSQCURVEFIT ITERATIONS
    options.MaxIter = options.MaxIter*4;
    options.MaxFunEvals= 600*4;
    param.kBoltz = 1.3807e-23; % Boltzmann's Constant [SI]
    param.c = 2.99792458e8; % speed of light [SI]

end
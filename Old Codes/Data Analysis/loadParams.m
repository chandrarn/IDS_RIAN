function[param, options] = loadParams(shot, line,useTree, s)
    x0shift = 0;%14;3; % shift x0 points %80
    y0shift = 0;-12; % shift y0 points %46, -28 
    
    if useTree % try using the tree, if that doesnt work, try importing callibration from file
        import MDSplus.*;
        if ~(shot>999999) % I ASSUME THAT ALL THE NODES I WANT EXIST ON HITSI AND HITSI3
            HitTree = Tree('hitsi',shot);
        else
            HitTree= Tree('hitsi3',shot);
        end


        Data = HitTree.getNode('\IDS_PEAKS');
        param.peaks = NATIVEvalue(Data.getData()); % [pixels] (n_channels x 5) chan number, x, y, sig_x, sig_y


        if ~isempty(x0shift)
            param.peaks(:, 2) = param.peaks(:, 2) + x0shift
        end

        param.peaks([11:15,42:43],2) = param.peaks([11:15,42:43],2) + 1; % ELIMINATE 21/1/16
        
        param.peaks(:,3) = param.peaks(:,3)+y0shift;

        
        Data = HitTree.getNode('\IDS_MASS');
        IonMass = (1.66e-27)*NATIVEvalue(Data.getData()); % [AMU] => [kg]
        param.IonMass = IonMass(line(1)); % select mass of ion of interest- for doublet, first one

        Data = HitTree.getNode('\CAL_LAMBDA');
        param.CalLam = NATIVEvalue(Data.getData());%-1e-9 * .15; % +1e-9 * 0.15[m]
        
        Data = HitTree.getNode('\IDS_LAMBDA');
        LineLam = NATIVEvalue(Data.getData()); % [m]
        param.LineLam = LineLam; % skip calibration lambda, which is first
        % for doublet, this automatically selects both lines
        % USed to be LineLam(line)
        
        Data = HitTree.getNode('\IDS_PIX_SP');
        param.PIX_SP = NATIVEvalue(Data.getData()); % [m] (1 x n_channels) wavelength spacing per channel
%         %param.PIX_SP = [1.58489019273755e-11;1.58355526981418e-11;1.58679347828107e-11;1.58534433956298e-11;1.58712475629411e-11;1.58596629278092e-11;1.58999964813235e-11;1.58779110472550e-11;1.58794286855083e-11;1.58777235776152e-11;1.58957301995632e-11;1.58868088277804e-11;1.59156736447004e-11;1.58835111788210e-11;1.59087665365898e-11;1.59142619817177e-11;1.58930307469100e-11;1.59161023403247e-11;1.59131704094126e-11;1.58755234582688e-11;1.59361285671668e-11;1.59005311017466e-11;1.59199443573188e-11;1.59252493486517e-11;1.59117630277443e-11;1.59074839844379e-11;1.58974315020601e-11;1.59460587315190e-11;1.58977718954227e-11;1.58965039192160e-11;1.59243696648995e-11;1.59022267147500e-11;1.59228599755694e-11;1.58749654494509e-11;1.58981881487834e-11;1.58942832113746e-11;1.59687404754893e-11;1.59242997100994e-11;1.59108180132415e-11;1.58702362365765e-11;1.58896167781856e-11;1.58646760261681e-11;1.59034969144890e-11;1.58567811849292e-11;1.58983180873064e-11;1.58976094108484e-11;1.58409122531724e-11;1.58756474070453e-11;1.59230899062054e-11;1.58958217966224e-11;1.58527035714640e-11;1.58996906569380e-11;1.58357751146298e-11;1.58463880580432e-11;1.58366689483885e-11;1.58405660687194e-11;1.58529948928585e-11;1.58362326815520e-11;1.58251773644243e-11;1.58521488651230e-11;1.58270727307333e-11;1.57815472879606e-11;1.57804243534196e-11;1.57746922948820e-11;1.57797395515402e-11;1.58207259984682e-11;1.57753817969537e-11;];
%         param.PIX_SP = [3.42249445827584e-12;1.14495371442319e-11;1.14683541661484e-11;1.14696251761017e-11;2.65664149860235e-12;1.14704880727463e-11;1.14663544135369e-11;1.14958991343689e-11;1.14830627501526e-11;1.14798467182889e-11;1.14815378467547e-11;1.14911275667289e-11;1.14863744933970e-11;1.15085319120920e-11;1.14880910884337e-11;1.14792149527205e-11;1.15083437381323e-11;1.14774173997467e-11;1.14913198799654e-11;1.15052654371721e-11;1.14992579860098e-11;1.14867066236053e-11;1.15047483621215e-11;1.14896594939180e-11;1.15205469873607e-11;1.15136501308894e-11;1.14921098383455e-11;1.14986123700954e-11;1.14849695732518e-11;1.15288069194312e-11;1.14852795515657e-11;1.14749453222570e-11;1.15012465201750e-11;1.15036048440263e-11;1.14866179857977e-11;1.14852140546411e-11;1.15253054840862e-11;1.14766240359109e-11;1.15313469665556e-11;1.14857672708946e-11;1.14909017946328e-11;1.15046635107152e-11;1.14736285500413e-11;1.14649002853401e-11;1.15263325548938e-11;1.14833058469721e-11;9.53039789151639e-12;1.15114949298141e-11;1.14627637123293e-11;1.15178848094988e-11;1.07053677818909e-11;5.05157844426785e-12;3.97644141114689e-12;5.88587090338135e-12;1.14500415629667e-11;1.14698732521001e-11;1.14679854366453e-11;1.14665412844822e-11;1.14514018399221e-11;1.14778776931451e-11;1.14477656461196e-11;1.15026649504550e-11;1.14212808742993e-11;1.14654138777293e-11;9.10379954814151e-12;9.14416824377639e-12;1.13849799567961e-11;]
        Data = HitTree.getNode('\IDS_REL_INT');
        param.REL_INT = NATIVEvalue(Data.getData()); % INVERSE relative intensity, 
                                       % ie: multiply by fit area to calibrate
        Data = HitTree.getNode('\IDS_IMPACTS');
        param.impacts=NATIVEvalue(Data.getData());
        mdsclose();
    else
        %cd('T:\IDS\Calibration\Calibration repository')
        temp = importdata(['T:\IDS\Calibration\Calibration repository\' int2str(shot) '.mat'] ); % load in params from file
        assignin('base','temp',temp);
        param.peaks = temp.PEAKS;
        param.IonMass = (1.66e-27)*temp.MASS(line(1));
        param.CalLam = temp.CAL_LAMBDA;
        param.LineLam = temp.LAMBDA(line);
        param.PIX_SP = temp.PIX_SP;
        param.REL_INT = temp.REL_INT;
        param.Impacts = temp.IMPACTS;
        

        
    end
    
    %param.peaks(9:27,2)=[41;44;47.5;50.5;54;57;61;64;67.5;71;74;78;81;88;91;95;99;102.5;106];
     %param.peaks(12,2)=89; param.peaks(15,2)=101.75; param.peaks(19,2)=118.5; param.peaks(8,2)=71.75;param.peaks(22,2)=136.5;  
    if length(line) == 2 % fit to a doublet
        param.dbl = 1; % flag to tell Gaussian fitting routine to fit doublet
    else
        param.dbl = 0; % flag indicating single line- not a doublet
    end

    param.shotRef = shot;
    size((param.CalLam - param.LineLam(line(end))));
    size(param.PIX_SP);
    % Adjust Center from calibration line to line of interest- NB: for DOUBLET,
    % center is arbitrarily chosen to be the smaller wavelength line ***
    % RNC: this only works if the calibration line and data line are in
    % almost the exact same place on the CCD. If care is not taken, this
    % will vastly over/under estimate the center location. Manually setting
    % this is generally easier.
    param.Center = param.peaks(:, 3);%-(param.CalLam - param.LineLam(line(end))) ./ param.PIX_SP %-10
    % HIGHLY TEMPORARY FIX
    
    % Establish Constants
    param.limits = [0 100; -40e3 40e3];
        % limits for sanitizing final output- throws out data  that exceeds the
        % above limits for [temp; vel]
    param.n_chan = size(param.peaks, 1);

    param.calcError = 1; % calculate uncertainty by LM method - presently broken
    param.ampThresh = 0; % [pixels], amplitude threshold TEMP SET TO ONE
    param.xTol = 0.5; % [pixels], tolerance of recentering on line in spatial dimension
    param.yTol = 3.5; % [pixels], tolerance for velocity offset
    param.xWing = 2; % real-space domain for fitting. 1 = 3 wide, 2 = 5 wide
    param.yWing = 7; % wavelength-space domain. 5 = 11 wide, etc...
    options = optimsetv61('lsqcurvefit'); % set options for curve fitting
    param.kBoltz = 1.3807e-23; % Boltzmann's Constant [SI]
    param.c = 2.99792458e8; % speed of light [SI]

end
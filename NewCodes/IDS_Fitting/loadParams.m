function[param, options] = loadParams(shot, line, hitsi3, useTree)
    
    % DELETE THIS IF YOU FIND IT:
    %shot = 150625006; % because we can't save to pdc3
    %spectraPix = 1.1660e-11;
    
    lineShift = [-14,0,11]; % attempt to manually fit to the lines;
    %y0 offset, set the cal line to OII
    y0 = 16;11;8;
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
        %param.peaks(find(param.peaks(:,1)==32):end,2) = param.peaks(find(param.peaks(:,1)==32):end,2) -1; % temp fix for lower array
%         param.peaks([11:18,20:21,27:28,42:49,52:55],2) = ... % X0 fix, manual
%             param.peaks([11:18,20:21,27:28,42:49,52:55],2) + 1; % ELIMINATE 21/1/16, for 151217026
%         param.peaks(:,2)=[10.4032211303711;18.8437366485596;22.8630275726318;27.0328025817871;31.2429351806641;35.3626976013184;39.3851165771484;43.4999961853027;47.8446655273438;51.8583030700684;57.1608161926270;61.2901458740234;65.2988357543945;69.4239120483398;73.5000000000000;77.8307952880859;81.9760971069336;86.3938293457031;89.6407775878906;94.8055343627930;99.4050903320313;103.091369628906;107.623367309570;112.115753173828;116.500030517578;120.619110107422;126.184967041016;130.282913208008;133.638137817383;168.182846069336;172.125213623047;176.224822998047;180.380523681641;184.577499389648;188.500000000000;192.791595458984;196.903701782227;200.937759399414;205.153228759766;209.132675170898;213.189437866211;218.525543212891;222.659332275391;226.430557250977;230.706726074219;234.934707641602;239.075592041016;243.093505859375;247.370925903320;250.500000000000;250.661163330078;250.856323242188;250.171112060547;250.353759765625;250.994598388672]' -2;
%         param.peaks(14:end,2) = param.peaks(14:end,2) +1.5;
%         param.peaks(20,2) = param.peaks(20,2)+.7;
%         param.peaks(21,2) = param.peaks(21,2)-.5;
%         param.peaks(23,2) = param.peaks(23,2)+.7;
        assignin('base','OrigPeaks',param.peaks);
        
        Data = Conn.get('\IDS_MASS');
        IonMass = (1.66e-27) * double(NATIVEvalue(Data.getFloatArray)); % [AMU] => [kg]
        param.IonMass = IonMass(line); % select mass of ion of interest- for doublet, first one

        Data = Conn.get('\CAL_LAMBDA');
        param.CalLam = double(NATIVEvalue(Data.getFloatArray)); % [m]
        Data = Conn.get('\IDS_LAMBDA');
        LineLam = double(NATIVEvalue(Data.getFloatArray)); % [m]
        param.LineLam = LineLam; % skip calibration lambda, which is first
        % for doublet, this automatically selects both lines
        % USed to be LineLam(line)
        
        Data = Conn.get('\IDS_PIX_SP');
        param.PIX_SP = double(NATIVEvalue(Data.getFloatArray)); % [m] (1 x n_channels) wavelength spacing per channel
        % TEMP FIX, DELETE IF FOUND
        %% HIGHLY TEMPORARY TESTING DoubletPixSp
        %param.PIX_SP = [1.53722578581739e-11,1.38212923051030e-11,1.19719624110047e-11,1.17026089394815e-11,1.14283499116489e-11,1.25750297840422e-11,1.15893539150753e-11,1.14511424653647e-11,1.13714463840290e-11,1.15120116997154e-11,1.11921237491581e-11,1.13650736053593e-11,1.14613876380423e-11,1.13651295468723e-11,1.13931223784949e-11,1.14011014334439e-11,1.13684065474829e-11,1.13937850860116e-11,1.13232692133044e-11,1.13878222899515e-11,1.13531509832043e-11,1.12861352857661e-11,1.13220242984979e-11,1.14159927485432e-11,1.14948175330188e-11,1.14250954586355e-11,1.20492604126363e-11,1.19362990840493e-11,1.35406221039772e-11,1.36145388892628e-11,1.55288012908145e-11,1.35180475722213e-11,1.25134368976061e-11,1.39718481814872e-11,1.39178413550214e-11,1.55286258367720e-11,1.28793465910828e-11,1.12072776816946e-11,1.14214709293663e-11,1.17558891254022e-11,1.13319878289780e-11,1.11769872104238e-11,1.11628302363138e-11,1.13269378956056e-11,1.13090290430226e-11,1.14515339837955e-11,1.12944663798623e-11,1.13623496038249e-11,1.12879767580595e-11,1.13671714665189e-11,1.13146427604977e-11,1.12778175360433e-11,1.13161226528156e-11,1.13478996347735e-11,1.18329562430502e-11;]';
        %param.PIX_SP = [6.10853001153053e-14,6.55634452560952e-12,1.12946931420782e-11,1.13256418815330e-11,1.12971473802400e-11,2.55079817058740e-13,1.04614939358760e-12,1.13137911626854e-11,1.13301424127941e-11,1.13296339497597e-11,1.13181224399386e-11,1.13262720393159e-11,1.13282715933229e-11,1.13199154503330e-11,1.13480970803359e-11,1.13279762097328e-11,1.13275644672364e-11,1.13328827528953e-11,1.13354676912299e-11,1.13501991621426e-11,1.13628624884301e-11,1.13318392619362e-11,1.13426876663099e-11,1.13561480935162e-11,1.13701541872217e-11,1.13147322048648e-11,6.16132445283354e-12,2.19847703761373e-12,4.00918411201151e-13,1.13563648313871e-11,1.13538887470013e-11,1.13480865593881e-11,1.13247470323049e-11,1.13287393316280e-11,1.13150645115702e-11,1.13065713786524e-11,1.13086664679395e-11,1.13080153872410e-11,1.13061652445182e-11,1.13217226876452e-11,1.13146427270259e-11,1.12992338080685e-11,1.12966144768122e-11,1.13193926339068e-11,1.13143703805354e-11,1.13047772550923e-11,1.13081308553339e-11,1.13174937731522e-11,1.12878355246104e-11,1.12993065073404e-11,1.12760079168827e-11,1.12773052474839e-11,1.12779408708449e-11,1.13016614893499e-11,1.12932149997692e-11]';
        %param.PIX_SP([1,2,6,7,27,28,29]) = mean(param.PIX_SP([3:5,9:26,30:end]));
        %param.PIX_SP = param.PIX_SP ;%* (spectraPix/mean(param.PIX_SP));
        % for 160517008
        %param.PIX_SP = param.PIX_SP(1:29);
        
        Data = Conn.get('\IDS_REL_INT');
        param.REL_INT = double(NATIVEvalue(Data.getFloatArray)); % INVERSE relative intensity, 
        %% ALSHO HIGHLY TEMPORARY 160517008
        %param.REL_INT = param.REL_INT(1:29);
        
                                       % ie: multiply by fit area to calibrate
        Data = Conn.get('\IDS_IMPACTS');
        param.impacts = double(NATIVEvalue(Data.getFloatArray));
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
        
        assignin('base','param',param)
        
    end
    
    %param.peaks(9:27,2)=[41;44;47.5;50.5;54;57;61;64;67.5;71;74;78;81;88;91;95;99;102.5;106];
     %param.peaks(12,2)=89; param.peaks(15,2)=101.75; param.peaks(19,2)=118.5; param.peaks(8,2)=71.75;param.peaks(22,2)=136.5;  
%     if length(line) == 2 % fit to a doublet
%         param.dbl = 1; % flag to tell Gaussian fitting routine to fit doublet
%     else
        param.dbl = 0; % flag indicating single line- not a doublet
%     end

    %offset the line
    %param.peaks(:,3) = param.peaks(:,3) +y0;
    % slightly more clevel offset, takes into account curvature. Kinda.
    %param.CalLam = 464.913e-9 -(y0 - mean(param.peaks(:,3)))*mean(param.PIX_SP)
    
    
    assignin('base','param',param);
%     x=imputdlg('TEST')
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
            param.Center(i, n) = param.peaks(i, 3)-(param.LineLam(3) - param.LineLam(line(n))) ./ temp; %mean(param.PIX_SP);%LineLam replaces CalLam, Mean(pixsp) underpredicts
         end
    end

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
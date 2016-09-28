% Aaron Hossack
% May 13th, 2013

%% UPDATED AARON HOSSACK, JAN 2015
% Making sure same code works on my linux machine. This code was updated
% with the best available calibration data which applies to all 53.5 and
% 68.5 kHz data. Original update was6 done for B.S. Victor's PoP paper.
%
%% UPDATED BY RIAN CHANDRA, MAY 2014
    % Code now inclues doublet calibration, can interface with new cine
    % file formats, and is updated to new MDSplus. SPR 2014.

% This is going to be an overarching calibration code that handles
% EVERYTHING.
%
% PHASE 1:
% Load in calibration movies for both upper and lower fibers, perform SVD,
% add first topos for both movies together, display.  User selects channel
% numbers and coordinates of first and last channel, and excludes dead
% channels.
%
% PHASE 2:
% Perform first fit for approximate coordinates of all channels.  Bin in
% WAVELENGTH direction and fit sum of all spatial direction Gaussians using
% the same fitting parameters for all.  This yields approximate center coordinates
% for all channels and estimate of spatial FWHM.
%
% Perform accurate, 2D Gaussian fit.  Loop through each channel, select
% three pixel columns nearest peak in real space and appropriate bounds in
% wavelength space.  Create mesh grid and fit to Gaussian.  Output 'peaks'
% and 'bounds' with option to save to tree.  Also outputs CENTER and FWHM.
%
%---------- Doublet Option
% PHASE 3:
% Step one of doublet fitting - bin in real direction to visually find
% separation of peaks
%
% PHASE 4:
% Optional doublet calibration - moves the specified distance away from the
% main line to a dimmer line and repeats the 2D Gaussian fit to obtain an
% extimate of PIX_SP.
%
%--------- Motor Option
% PHASE 5:
% Motor Calibration for PIX_SP - based on slant of earlier calibration
% lines, calculates function for reestimating center at any wavelength
% position on CCD
%
% PHASE5B:
% The grating motor is very linear and consistent, but apparently
% SIGNIFICANTLY not close to its claimed speed.  This section corrects for
% this by loading in a plasma shot and fitting to two or more lines, then
% adjusting all the PIX_SP values.
%
%---------
%
% PHASE 6:
% Find FWHM and plot FWHM, PIX_SP, REL_INT
%
% PHASE 7:
% Discard marginal or bad channels
%
% PHASE 8:
% Save selected settings to tree.
%
clear all; close all; clc;
addpath('~/IDS/Matlab/');
try
    addAllThePaths;
end

%% Settings

% SAVING FIGURES --------
% PHASE 1
approxFit.save = 0;
approxFit.file = 'T:\RChandra\A-A-Ron Code\test';

% PHASE 2
indivFit.save = 0;
indivFit.chan = 24;
indivFit.file = ['T:\RChandra\A-A-Ron Code\example_cal_fit_' ...
    num2str(indivFit.chan)];

allFit.save = 0;
allFit.file = 'T:\RChandra\A-A-Ron Code\\test';

% PHASE 3
dblBinFit.save = 0;
dblBinFit.file = [];

%%
% GENERAL -----------------------------------------------------------------
xWing = 1; % real space domain (1 => 3, 2 => 5)

% PHASE 1 -----------------------------------------------------------------
doPHASE1 = 1;
shot1 = 13041901; % always use this one
%shot2 = 0; % set 'shot2' to zero if only using one fiber
% shot1 = 13041901; % always use this one
shot2 = 13041903; % set 'shot2' to zero if only using one fiber
xlim = [0 256]; % pixels for initial plotting, leave empty if unsure
ylim = [48 68]; % pixels for initial plotting, leave empty if unsure

% PHASE 2 -----------------------------------------------------------------
doPHASE2 = 1;

chanNums = [1:4, 6:22, 24:36, 37:72];%, 72:86]; % 
firstCenter = [12, 58]; % Center position of first channel, [real, wavelength]
lastCenter = [254, 58]; % Center position of last channel, [real, wavelength]
brightWing = 4; % number of pixels in wavelength space for Gaussian fitting domain
force = [12, 77, 83, 193, 229]; % force finding channel(s) at specific x location(s)
remove = [84]; % pixel number in channel space to ignore a local max

% ------- 14.5 kHz, May/June 2013
% % chanNums = [1:4, 6:21, 23:72]; % 5 and 22 are dead
% chanNums = [1:4, 6:72]; % 5 dead, 22 weak but including - 5/31/13
% firstCenter = [37, 55]; % Center position of first channel, [real, wavelength]
% lastCenter = [348, 55]; % Center position of last channel, [real, wavelength]
% brightWing = 5; % number of pixels in wavelength space for Gaussian fitting domain
% force = [];

% PHASE 3 -----------------------------------------------------------------
doPHASE3 = 1;
DoubletShot1 = 13030401;
DoubletShot2 = 13030402;


nPix2Bin = 20; % number of pixels in real space (to the right of channel 1) to bin for display

% PHASE 4 -----------------------------------------------------------------
doPHASE4 = 1;
largerLam = 463.723e-9;
% smallerLam = 460.957e-9;
smallerLam = 460.463e-9; % This is a kludgy correction to get PIX-SP to come out to the correct value - ACH, 2-3-15

% PHASE 5 -----------------------------------------------------------------
% MUST MANUALLY CORRECT TIME BASE !
% BOTH MOVIES MUST HAVE SAME TIME INTERVAL !
doPHASE5 = 0;
shot3 = 14062504; % Motor Calibration Shot
shot4 = 14062505; % Optional second motor calibration shot
motorSpeed = 1; % [nm per second]

% PHASE 5 B ---------------------------------------------------------------
doPHASE5B = 0;
shotPlas = 140626021; % Plasma shot to analyze
channel = 21; % Channel Number
timePt = 50; % time point
pixelNums = [39 55 65]; % pixel numbers of peaks to within 1 pixel
plasmaLams = 1e-9 * [465.025, 464.913, 464.742]; % corresponding line wavelengths
calLam = 1e-9 * 465.025; % [m]
yWing = 7; % similar to 'brightWing', but allowing different value for plasma lines
factor = 0.7; % estimated correction to motor speed

% PHASE 6 -----------------------------------------------------------------
doPHASE6 = 1;

% PHASE 7 -----------------------------------------------------------------
doPHASE7 = 1;
discard = [22, 53]; % added 22 on 4-27-15, ACH

% ------- 14.5 kHz, May/June 2013
% discard = [22, 53, 54];

% PHASE 8 -----------------------------------------------------------------
doPHASE8 = 1; % VERY IMPORTANT
hitsi3 = 0;

%% NB !!!! Different CAL_LAMBDA for 53.5 and 68.5 kHz !!!
% saveToShots =[129225:129226]; % adding a couple more 68.5 kHz shots, 4-22-15, ACH

saveToShots = [128457:128597, 129213:129215, 129225:129226]; 
                    % 53.5 kHz Toridal fiber, 71 degree re-enterant port, polidal fiber on 45 degree axial port, x-side
                    % ALSO valid for 68.5 kHz
impactsFile = '/home/aaron/IDS/Geometry/impacts1.mat'; 

% saveToShots = [129510:129534, 129544:129571]; % mohawk port, orthogonal to midplane, and axial port. 'impacts4'

%saveToShots = [129438:129453, 129460:129480, 129485:129500, 129577:129600, ...
%    129611:129627, 129642:129661, 129662:129685, 129688:129711, ...
%    129719:129746, 129747:129759, 129800:129802]; % mohawk port, chords in midplane, axial port, 'impacts2'

% saveToShots = [129787:129799]; % axial ports facing each other

% saveToShots = [129807:129824]; % reentrant port at 71 degrees and axial port, 'impacts1'

stt.CAL_LAMBDA = 1;
stt.LAMBDA = 1;
stt.VOLTAGE = 1;
stt.MASS = 1;
stt.PEAKS = 1;
stt.REL_INT = 1;
stt.PIX_SP = 1;
stt.IMPACTS = 1;

CAL_LAMBDA = 1e-9 * 464.913; % [m], should be effective wavelength that matches PEAKS centers 68.5 kHz !!!
% CAL_LAMBDA = 1e-9 * 465.1096; % [m] % 53.5 kHz !!!

LAMBDA = 1e-9 * [465.025, 464.913, 464.742]; % [m], plasma lines in descending wavelength
VOLTAGE = 900; % image intensifier voltage, '0' if not using
MASS = [12, 16, 12]; % [AMU], corresponding to LAMBDA

% 'impacts1' = long fiber in old reentrant port, toroidal midplane diameter,
% short fiber in axial port at 45 degrees, X side (June 2012 - May 2013)
% 'impacts2' = long fiber in mohawk reentrant port, fan in midplane, axial fiber 
% on poloidal section.

%% PHASE 1 - Load Data, SVD, Display

if doPHASE1
    [data, X, Y, n_pix, n_chan] = calSVD(shot1, shot2);

    %% Plot First Topo

    S = get(0,'ScreenSize');
    fntsz = 24;

    h1 = figure('Visible','on','Name','Calibration Data','Position',...
        [S(3)/12, S(4)/8, 5*S(3)/6 3*S(4)/4], 'Color', [1 1 1]);
    ax(1) = axes('Parent', h1, 'Position', [.1 .11 .8 .33], 'FontSize', fntsz);
    h3 = surf(X, Y, data);
    hold on;
    shading interp;
    colormap jet;
    colorbar;
    grid on;
    view([0 90]);
    set(gca, 'FontSize', fntsz);
    if length(xlim) == 2
        set(gca, 'XLim', xlim);
    else
        set(gca, 'XLim', [1 n_chan]);
    end
    if length(ylim) == 2
        set(gca, 'YLim', ylim);
    else
        set(gca, 'YLim', [1 n_pix]);
    end

    title('BD Filtered Calibration Data');
    xlabel('Pixel Number (real space)');
    ylabel({'Pixel Number'; '(wavelength space)'});
end

%% PHASE 2 - Approximate fit to find all channel centers

if doPHASE2
    
    % Find approximate center positions of all channels

    peaks = fitAllChans(data, chanNums, firstCenter, lastCenter, brightWing, approxFit, force, remove);

    % Fit each peak to 2D Gaussian individually

    [PEAKS, REL_INT, par, fits] = calGauss2D(peaks, data, brightWing, xWing, indivFit);

    % Plot Fits ------------------------------------

    figure(h1);
    ax(2) = axes('Parent', h1, 'Position', [.1 .57 .8 .33], 'FontSize', fntsz);
    for n = 1:size(peaks, 1) % loop over channels
        h3 = surf(fits.Xf(:,:,n), fits.Yf(:,:,n), fits.Zf(:,:,n));
        hold on;
    end
    shading interp;
    colormap jet;
    colorbar;
    view([0 90]);
    set(gca, 'FontSize', fntsz);
    set(gca, 'XTick', []);
    linkaxes(ax, 'xy');

    title('Reconstructed Fits');
%     xlabel('Pixel Number (real space)');
    ylabel({'Pixel Number'; '(wavelength space)'});
    pause(0.1);
    gcf;

    % SAVE
    if allFit.save
        fig_save = getframe(h1);
        [Xfig, mapfig] = frame2im(fig_save);
        imwrite(Xfig, [allFit.file '.png']);
    end

end

%% PHASE 3 - Doublet Calibration - bin in real space to estimate offset of second peak

if doPHASE3
    % INCOMPLETE
    [data, X, Y, n_pix, n_chan] = calSVD(DoubletShot1, DoubletShot2);%Shot3&4 probably.
    %We can plot this by copying all of the code over from the plot BD
    %section from Phase 1.  May want to make it a seperate function and
    %just call it twice.
    
    
    ypeaks = fitAllChansWavelength(data);
    
    % find y bounds for plotting
    ylim2 = [];
    xlim2 = [];
    ylim2 = [(ypeaks(1, 2) - 15) (ypeaks(2, 2) + 15)]; % in theory, this will give good bounds.
    
    S = get(0, 'ScreenSize');
    fntsz = 24;

    h1 = figure('Visible', 'on', 'Name', 'Calibration Data2', 'Position',...
        [S(3)/12, S(4)/8, 5*S(3)/6 3*S(4)/4], 'Color', [1 1 1]);
    %ax(1) = axes('Parent', h1, 'Position', [.1 .11 .8 .33], 'FontSize', fntsz);
    ax(1) = axes('Parent', h1, 'FontSize', fntsz);
    h3 = surf(X, Y, data);
    h3 = surf(X, Y, data);
    hold on;
    plot3(0:250, ypeaks(1, 2)*ones(1, 251), .025*ones(1, 251), 'linewidth', 3, 'Color', 'black');
    plot3(0:250, ypeaks(2, 2)*ones(1, 251), .025*ones(1, 251), 'linewidth', 3, 'Color', 'black');
    shading interp;
    colormap jet;
    colorbar;
    grid on;
    view([0 90]);
    set(gca, 'FontSize', fntsz);
    if length(xlim2) == 2
        set(gca, 'XLim', xlim2);
    else
        set(gca, 'XLim', [1 260]);
    end
    if length(ylim2) == 2
        set(gca, 'YLim', ylim2);
    else
        set(gca, 'YLim', [1 720]);
    end

    title('BD Filtered Doublet Calibration Data');
    xlabel('Pixel Number (real space)');
    ylabel({'Pixel Number'; '(wavelength space)'});
end

%% PHASE 4 - Doublet Calibration - Calculate PIX_SP from doublet fitting

if doPHASE4
    % INCOMPLETE
    %Test doing rough channels fit for x guess values.
    [m, n] = size(data);
%     BottomXpeaks = fitAllChans(data(1:round(m/2),:), chanNums, firstCenter, lastCenter, brightWing, approxFit, force);
%     TopXpeaks = fitAllChans(data(round(m/2):end,:), chanNums, firstCenter, lastCenter, brightWing, approxFit, force);
    %take ypeaks, use that as the y-guess for 2dGauss, use peaks as xguess
    peaks(:, 1) = chanNums;
    peaks(:, 2) = par(:,2) + 1;%take the 2dgauss x-values from phase 2 
    peaks(:, 3) = ypeaks(1, 2);%do bottom line first (ypeaks(1) is lower index)
    [~, ~, par1, fits1] = calGauss2D(peaks, data, brightWing, xWing, indivFit);
    peaks(:, 3) = ypeaks(2,2);
    %peaks(:,2)=par1(:,2); Tested: didnt improve accuracy much.
    [~, ~, par2, fits2] = calGauss2D(peaks, data, brightWing, xWing, indivFit);
    %Test rough fit x
%     peaks(:,2)=BottomXpeaks(:,2);
%     peaks(:,3)=ypeaks(1,2);
%     [PEAKS, REL_INT, par3, fits1] = calGauss2D(peaks, data, brightWing, xWing, indivFit);
%     peaks(:,2)=TopXpeaks(:,2);
%     peaks(:,3)=ypeaks(2,2);
%     [PEAKS, REL_INT, par4, fits1] = calGauss2D(peaks, data, brightWing, xWing, indivFit);
    
    WavelengthOffset = par2(:, 3) - par1(:, 3);
    %Gauss fit using singlet gauss x values for x guess
    scatter3(par2(:, 2), par2(:, 3), 0.03 * ones(1, size(par2(:, 3), 1)), 'r', 'Marker', '*');
    scatter3(par1(:, 2), par1(:, 3), 0.03 * ones(1, size(par1(:, 3), 1)), 'r', 'Marker', '*');
    %Guess value, singlet gauss x and rough y
    scatter3(par(:, 2) + 1, ypeaks(1, 2) * ones(1, size(par1(:, 3), 1)), 0.03 * ones(1, size(par1(:, 3), 1)), 'g', 'Marker', '*');
    scatter3(par(:, 2) + 1, ypeaks(2, 2) * ones(1, size(par1(:, 3), 1)), 0.03 * ones(1, size(par1(:, 3), 1)), 'g', 'Marker', '*');
    %rough x and rough y
%     scatter3(BottomXpeaks(:,2),ypeaks(1,2)*ones(1,size(par1(:,3))),.03*ones(1,size(par1(:,3))),'y','Marker','*');
%     scatter3(TopXpeaks(:,2),ypeaks(2,2)*ones(1,size(par1(:,3))),.03*ones(1,size(par1(:,3))),'y','Marker','*');
    %gauss fit with rough X & y
%     scatter3(par3(:,2),par3(:,3),.03*ones(1,size(par3(:,3))),'c','Marker','*');
%     scatter3(par4(:,2),par4(:,3),.03*ones(1,size(par4(:,3))),'c','Marker','*');
    iterator = 1; %create two vectors of points, to plot line segments
    Xvector = 0;
    Yvector = 0;
    for points=1:size(peaks(:, 2), 1)
        Xvector(iterator) = par1(points, 2);
        Xvector(iterator + 1) = par2(points, 2);
        Xvector(iterator + 3) = 0;
        Yvector(iterator) = par1(points, 3);
        Yvector(iterator + 1) = par2(points, 3);
        Yvector(iterator + 3) = 0;
        iterator = iterator + 3;
    end
    Xvector(3:3:(size(Xvector, 2))) = NaN;
    Yvector(3:3:(size(Xvector, 2))) = NaN;
    plot3(Xvector, Yvector, 0.03*ones(1, size(Yvector, 2)), 'color', 'red');
    
    
    %%% Divide wavelength by pixel offset
    DifferenceInWavelength = ((largerLam) - (smallerLam)) .* ones(size(WavelengthOffset));
    PIX_SP = DifferenceInWavelength ./ WavelengthOffset;
    
%     % Plot Gauss Fits ------------------------------------
% 
%     h2 = figure('Visible','on','Name','Gauss doublet fits','Position',...
%         [S(3)/12, S(4)/8, 5*S(3)/6 3*S(4)/4], 'Color', [1 1 1]);
%     ax(1) = axes('Parent', h2, 'Position', [.1 .57 .8 .33], 'FontSize', fntsz);
%     for n = 1:size(peaks, 1) % loop over channels
%         h2 = surf(fits2.Xf(:,:,n), fits2.Yf(:,:,n), fits2.Zf(:,:,n));
%         hold on;
%     end
%     shading interp;
%     colormap jet;
%     colorbar;
%     view([0 90]);
%     set(gca, 'FontSize', fntsz);
%     set(gca, 'XTick', []);
%     linkaxes(ax, 'xy');
% 
%     title('Top Line ');
% %     xlabel('Pixel Number (real space)');
%     ylabel({'Pixel Number'; '(wavelength space)'});
%     pause(0.1);
%     %axis([0,256,ypeaks(1,2)-5,ypeaks(1,2)+5]);
%     set(gca,'YLim',[(ypeaks(1,2)-5) (ypeaks(1,2)+5)]);
%     gcf;
%     
%     
%     
%     %plot second set of figures
%     figure(h2);
%     ax(2) = axes('Parent', h1, 'Position', [.1 .57 .8 .33], 'FontSize', fntsz);
%     for n = 1:size(peaks, 1) % loop over channels
%         h3 = surf(fits1.Xf(:,:,n), fits1.Yf(:,:,n), fits1.Zf(:,:,n));
%         hold on;
%     end
%     shading interp;
%     colormap jet;
%     colorbar;
%     view([0 90]);
%     set(gca, 'FontSize', fntsz);
%     set(gca, 'XTick', []);
%     linkaxes(ax, 'xy');
% 
%     title('Bottom line');
% %     xlabel('Pixel Number (real space)');
%     ylabel({'Pixel Number'; '(wavelength space)'});
%     pause(0.1);
%     gcf;

end

%% PHASE 5 - Motor Calibration for PIX_SP

if doPHASE5
    
    % Step 1: Find parameters to reestimate center anywhere in wavelength
    % domain of CCD
    
%     slope = findSlant(PEAKS); % Sidelined for now- assume fibers aligned
%     with CCD
    
    % Step 2: Loop through time and fit every channel
    
    PIX_SP = fitMotor(shot3, shot4, PEAKS, motorSpeed, brightWing, xWing);
    
end

%% PHASE 5 B - Correct Motor Calibration with Plasma Lines

if doPHASE5B
    %input('HALTING PHASE FIVE B')
    PIX_SP = calPlasma(shotPlas, PIX_SP, PEAKS, channel, timePt, ...
        pixelNums, plasmaLams, calLam, xWing, yWing, factor);

end

%% PHASE 6 - Convert Fitting Parameters, plot FWHM, PIX_SP, and REL_INT

if doPHASE6
    
    FWHM = 2.35482 * PIX_SP .* PEAKS(:, 5); % convert sigma to FWHM
    
    h1 = figure('Visible','on','Name','Approximate Fit','Position',...
        [S(3)/12, S(4)/6, 5*S(3)/6 2*S(4)/3], 'Color', [1 1 1]);
    
    % REL_INT
    h2 = axes('Parent', h1, 'Position', [.1 .1 .8 .25], 'FontSize', fntsz);
    h3 = plot(PEAKS(:, 1), REL_INT, 'ob');
    hold on;
    grid on;
    set(gca, 'XLim', [PEAKS(1, 1) PEAKS(end, 1)]);
    xlabel('Channel Number');
    ylabel('REL INT');
    
    % FWHM
    h2 = axes('Parent', h1, 'Position', [.1 .4 .8 .25], 'FontSize', fntsz);
    h3 = plot(PEAKS(:, 1), FWHM, 'ob');
    hold on;
    grid on;
    set(gca, 'XLim', [PEAKS(1, 1) PEAKS(end, 1)]);
    set(gca, 'XTickLab', []);
    ylabel('FWHM [m]');
    
    % PIX_SP
    h2 = axes('Parent', h1, 'Position', [.1 .7 .8 .25], 'FontSize', fntsz);
    h3 = plot(PEAKS(:, 1), PIX_SP, 'ob');
    hold on;
    grid on;
    set(gca, 'XLim', [PEAKS(1, 1) PEAKS(end, 1)]);
    set(gca, 'XTickLab', []);
    ylabel('PIX SP [m]');
    
end

%% PHASE 7 - Discard Bad Channels

if and(doPHASE7, length(discard >= 1))

    bad = [];
    for n = 1:length(discard)
        b = find(PEAKS(:, 1) == discard(n)); % index of nth bad channel
        bad = [bad b]; % all indices of bad channels
    end

    c = PEAKS(:, 1); % alternate array of channel numbers
    c(bad) = zeros(1, length(discard)); % zero out bad channels so I can use 'find'
    ind = find(c); % finds all non-zero elements of c

    % Use 'ind' to trim arrays

    PEAKS = PEAKS(ind, :);
    REL_INT = REL_INT(ind);
    PIX_SP = PIX_SP(ind);
end

%% Retrieve 'impacts', print all data before saving

try load(impactsFile); end;
IMPACTS = impacts(PEAKS(:, 1)); % select only channel numbers being used

try
    REL_INT;
catch % not running the whole code
    REL_INT = NaN;
    PEAKS = NaN;
    PIX_SP = NaN;
    stt.REL_INT = 0; % make damn sure I don't write NaNs to the tree
    stt.PEAKS = 0;
    stt.PIX_SP = 0;
end

saveToShots
CAL_LAMBDA
LAMBDA
VOLTAGE
MASS
PEAKS
REL_INT
IMPACTS
PIX_SP



%% PHASE 8 - Save Settings to Tree

if doPHASE8
    try
        saveToTree3(saveToShots, stt, CAL_LAMBDA, LAMBDA, VOLTAGE, MASS, PEAKS, REL_INT, PIX_SP, IMPACTS, hitsi3);
    
    catch 
        cd('T:\IDS\Calibration\Calibration repository')
        fsave.CAL_LAMBDA = CAL_LAMBDA;
        fsave.LAMBDA = LAMBDA;
        fsave.VOLTAGE = VOLTAGE;
        fsave.MASS = MASS;
        fsave.PEAKS = PEAKS;
        fsave.REL_INT = REL_INT;
        fsave.PIX_SP = PIX_SP;
        fsave.IMPACTS = IMPACTS;

        for i = saveToShots
            save(num2str(i),'fsave')
        end
    end
end

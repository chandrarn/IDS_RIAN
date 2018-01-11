% Absolute error correction
% Predicts gaussian position based on initial fit
% Based on calMain Phases 1,2,5a
% Calculates deviation from modeled frame. Store point pairs, Int and Err,
% for each pixel


%function [Int,Err] = absError(motorCalShot,trimTime,binChanMotor,xWing,oldFileType,param)
% Predefine variables
motorCalShot = 13053104;
trimTime = [1:22]';
binChanMotor = 166;
xWing = 2;
oldFileType=1;
channel = 30;
param = param;
% Phase 1,2
shot1 = 13053101;
shot2 = 13053102;
chanNums = [1:4, 6:72]; % Nothing wrong with 54, if lens is removed durring calibration
firstCenter = [4,55]; % Center position of first channel, [real, wavelength]
lastCenter = [309,56]; % Center position of last channel, [real, wavelength]
brightWing = 5; % number of pixels in wavelength space for Gaussian fitting domain
force = []; % force finding channel(s) at specific x location(s)
remove = [];
% Phase 5
shot3 = 170504008; % Motor Calibration Shot
shot4 = 170504006; % Optional second motor calibration shot
motorSpeed = 0.0509; % [nm per second]
trim5 = [1:22]; % Trim calibration movie to just bright line

%% Load Data
% calMain Phases 1-2
[data, X, Y, n_pix, n_chan] = calSVD(shot1, shot2);
PEAKS = fitAllChans(data, chanNums, firstCenter, lastCenter, brightWing, [], force, remove);
% Fit each peak to 2D Gaussian individually
[PEAKS, ~, ~, ~] = calGauss2D(PEAKS, data, brightWing, xWing, []);
%% Fit Gaussians
% calMain Phase 5
[~,par,data] = fitMotor(shot3, shot4, PEAKS, motorSpeed, brightWing, xWing,trim5);

%% Calculate slope
% calculate slope of all parmeters
parSlope = zeros(n_chan,6,2);
for n = 1:n_chan
    for i = 1:6
        [parSlope(n,i,:),~] = polyfit(trim5,par(:,n,i),1); % fit param vs timepoint
    end
end

% Visualize the parameters
figure;
for i = 1:6
subplot(2,3,i);
plot(trim5,par(:,:,i),'*');
hold on; plot(trim5.*mean(parSlope(:,i,1)) + mean(parSlope(:,i,2)),'-r');
title(['(' num2str(mean(parSlope(:,i,1))) '\pm' num2str(std(parSlope(:,i,1)))...
    ')x + ' '(' num2str(mean(parSlope(:,i,2))) '\pm' num2str(std(parSlope(:,i,2))) ')']);
end

 
 %% Apply static calibration
 [X, Y] = meshgrid(1:size(data,2), 1:size(data,3));
 x(:, 1) = X(:);
x(:, 2) = Y(:);
model = zeros(1,5);
Int = zeros(size(data));
Err = zeros(size(data));
 for t = trim5 % loop over timePoints
     % Model frame for time t
      modelFrame = zeros(size(data));
     for n = 1:param.n_chan
         model(:) = [parSlope(n,:,1).*t + parSlope(n,:,2)];
         modelFrame = modelFrame + singletGauss2D(model, x);
     end
     Int(t,:,:)=modelFrame;
     Err(t,:,:) = data(t,:,:)-Int;
 end
 
 
%end
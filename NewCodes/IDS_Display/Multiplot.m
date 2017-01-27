% Multiplot
% Note: need to manually save the figure as .eps, otherwise it saves in
% black and white. You know, because fuck all.
%close all; 
clear all;% clc;
addpath('~/IDS/Matlab/');
addpath('T:\RChandra\Sine_fit\');
addpath('T:\IDS\General Matlab\');
%addAllThePaths;
lines = {'O II', 'C III', 'O II','C III'};

%% Input Settings

%% IDS DATA
%% 0-120-240
% High performance
% % in(1).shot = 160728011;150625998;
% % in(1).line =2; % line # NB: 1 is C III, 2 is O II, 3 is C III !
% % in(1).legend = [num2str(in(1).shot) ' ' lines{in(1).line}];
% % in(1).legend = [num2str(in(1).shot) ': +65.5kA'];
% % in(1).color = {'r';'r'};%[225,105,0]./255};
% % in(1).style = {'-','--'};
% % in(1).error = 0; % 1 / 0 for errorbars
% % in(1).velShift = -5; % SHIFT VELOCITY
% % in(1).intScale = 2; % scale factor for intensity
% % in(1).timeShift = 0; % ms, shift time base
% % in(1).timeScale = 1e-3; % scale timebase to put into ms
% % in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
% % in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
% % in(1).doubleplot = [1];[1:23; 24,26:47]; % plot coorespoinding impacts
% % in(1).fftPlot = [1]; % FFT of signal, n frequencies
% % in(1).AnalysisTitle=[ 'HIT-SI3: 0-120-240 Phasing, Gain 2.9, ' lines{in(1).line}];
% % in(1).phaseShift = 2*pi -pi/2;
% % 160728
% %in(2)=in(1);
% in(1).shot = 160728013;
% in(1).line=2;
% in(1).color = {'r';'r'};%[66, 188, 244]./255};
% in(1).style = {'-','--'};
% in(1).legend = [num2str(in(1).shot) ': +65.6kA'];
% in(1).phaseShift = -pi/2;
% in(1).error = 0; % 1 / 0 for errorbars
% in(1).velShift = -5; % SHIFT VELOCITY
% in(1).intScale = 2; % scale factor for intensity
% in(1).timeShift = 0; % ms, shift time base
% in(1).timeScale = 1e-3; % scale timebase to put into ms
% in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
% in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(1).doubleplot = [1];[1:23; 24,26:47]; % plot coorespoinding impacts
% in(1).fftPlot = []; % FFT of signal, n frequencies
% in(1).AnalysisTitle=[ 'HIT-SI3: 0-120-240 Phasing, Gain 2.9, ' lines{in(1).line}];
% in(1).phaseShift = 2*pi -pi/2 ;% The +pi is to overlay with the positive shot
% in(2)=in(1);
% in(2).phaseShift = 2*pi -pi/2-pi;
% in(2).shot = 160728012;
% in(2).color = {'b';'b'};%[182, 244, 66]./255};
% in(2).legend = [num2str(in(2).shot) ': -52.3kA'];
% %Supression of Dissenting data:
% phaseSupress(1,:,1)=1;
% phaseSupress(1:2,:,2)=0;
% flowSupress(1,:)=[1:3];
% flowSupress(2,:)=[1:3];
% dispSupress(1,:,2)=[1:5];
% dispSupress(2,:,2)=[1:4,4];
% tempSupress(1,:,1)=[2,2,2];
% tempSupress(1,:,2)=[1:2,5];
% tempSupress(2,:,2)=[1,3,5];


% Low Performance
% in(1).shot = 160525017;150625998;
% in(1).line =2; % line # NB: 1 is C III, 2 is O II, 3 is C III !
% in(1).legend = [num2str(in(1).shot) ' ' lines{in(1).line}];
% in(1).legend = [num2str(in(1).shot) ': -41.9kA'];
% in(1).color = {'r';'r'};%[225,105,0]./255};
% in(1).style = {'-','--'};
% in(1).error = 0; % 1 / 0 for errorbars
% in(1).velShift = -5; % SHIFT VELOCITY
% in(1).intScale = 2; % scale factor for intensity
% in(1).timeShift = 0; % ms, shift time base
% in(1).timeScale = 1e-3; % scale timebase to put into ms
% in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
% in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(1).doubleplot = [1];[1:23; 24,26:47]; % plot coorespoinding impacts
% in(1).fftPlot = [1]; % FFT of signal, n frequencies
% in(1).AnalysisTitle=[ 'HIT-SI3: 0-120-240 Phasing, ' lines{in(1).line}];
% in(1).phaseShift = -pi/2 + pi;
% % % 160728
% in(2)=in(1);
% in(2).shot = 160525018;
% in(2).line=2;
% in(2).color = {'b';'b'};%[66, 188, 244]./255};
% in(2).legend = [num2str(in(2).shot) ': -45.5kA'];
% in(2).phaseShift = -pi/2+pi;
% in(3)=in(1);
% in(3).shot = 160525019;
% in(3).color = {'g';'g'};%[182, 244, 66]./255};
% in(3).legend = [num2str(in(3).shot) ': -43.4kA'];
% in(3).phaseShift = -pi/2+pi;


% in(1).shot = 160622017;150625998;
% in(1).line =2; % line # NB: 1 is C III, 2 is O II, 3 is C III !
% in(1).legend = [num2str(in(1).shot) ' ' lines{in(1).line}];
% in(1).legend = [num2str(in(1).shot) ': -47.1kA'];
% in(1).color = {'r';'r'};%[225,105,0]./255};
% in(1).style = {'-','--'};
% in(1).error = 0; % 1 / 0 for errorbars
% in(1).velShift = -5; % SHIFT VELOCITY
% in(1).intScale = 2; % scale factor for intensity
% in(1).timeShift = 0; % ms, shift time base
% in(1).timeScale = 1e-3; % scale timebase to put into ms
% in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
% in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(1).doubleplot = [1];[1:23; 24,26:47]; % plot coorespoinding impacts
% in(1).fftPlot = [1]; % FFT of signal, n frequencies
% in(1).AnalysisTitle=[ 'HIT-SI3: 0-120-240 Phasing, ' lines{in(1).line}];
% in(1).phaseShift = 2*pi-pi/2;
% % % 160728
% in(2)=in(1);
% in(2).shot = 160622018;
% in(2).line=2;
% in(2).color = {'b';'b'};%[66, 188, 244]./255};
% in(2).legend = [num2str(in(2).shot) ': -47.1kA'];
% in(2).phaseShift = -pi/2;
% in(3)=in(1);
% in(3).shot = 160622019;
% in(3).color = {'g';'g'};%[182, 244, 66]./255};
% in(3).legend = [num2str(in(3).shot) ': -50.3kA'];
% in(3).phaseShift = -2*pi-pi/2;



%% 0-60-120
% in(1).shot = 160728021;150625998;
% in(1).line =1; % line # NB: 1 is C III, 2 is O II, 3 is C III !
% in(1).legend = [num2str(in(1).shot) ' ' lines{in(1).line}];
% in(1).legend = [num2str(in(1).shot) ': -48.3kA'];
% in(1).color = {'r';'r'};%[225,105,0]./255};
% in(1).style = {'-','--'};
% in(1).error = 0; % 1 / 0 for errorbars
% in(1).velShift = -5; % SHIFT VELOCITY
% in(1).phaseShift = -pi/2; % Shift Phase
% in(1).intScale = 2; % scale factor for intensity
% in(1).timeShift = 0; % ms, shift time base
% in(1).timeScale = 1e-3; % scale timebase to put into ms
% in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
% in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(1).doubleplot = [1];[1:23; 24,26:47]; % plot coorespoinding impacts
% in(1).fftPlot = [1]; % FFT of signal, n frequencies
% in(1).AnalysisTitle=[ 'HIT-SI3: 0-60-120 Phasing, ' lines{in(1).line}];
% % 160728
% in(2)=in(1);
% in(2).shot = 160728023;
% in(2).line=1;
% in(2).color = {'b';'b'};%[66, 188, 244]./255};
% in(2).legend = [num2str(in(2).shot) ': -50.8kA'];
% in(3)=in(1);
% in(3).shot = 160728024;
% in(3).color = {'g';'g'};%[182, 244, 66]./255};
% in(3).legend = [num2str(in(3).shot) ': +55.2kA'];
% in(3).phaseShift = 2*pi -pi/2;

% Low Performance
%  in(1).shot = 160615024;150625998;
% in(1).line =2; % line # NB: 1 is C III, 2 is O II, 3 is C III !
% in(1).legend = [num2str(in(1).shot) ' ' lines{in(1).line}];
% in(1).legend = [num2str(in(1).shot) ': -47.1kA (+)'];
% in(1).color = {'r';'r'};%[225,105,0]./255};
% in(1).style = {'-','--'};
% in(1).error = 0; % 1 / 0 for errorbars
% in(1).velShift = -5; % SHIFT VELOCITY
% in(1).intScale = 2; % scale factor for intensity
% in(1).timeShift = 0; % ms, shift time base
% in(1).timeScale = 1e-3; % scale timebase to put into ms
% in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
% in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(1).doubleplot = [1];[1:23; 24,26:47]; % plot coorespoinding impacts
% in(1).fftPlot = [1]; % FFT of signal, n frequencies
% in(1).AnalysisTitle=[ 'HIT-SI3: 0-60-120 Phasing, ' lines{in(1).line}];
% in(1).phaseShift = -pi/2;
% % 160728
% in(2)=in(1);
% in(2).shot = 160615025;
% in(2).line=2;
% in(2).color = {'b';'b'};%[66, 188, 244]./255};
% in(2).legend = [num2str(in(2).shot) ': -47.1kA (-)'];
% in(2).phaseShift = 2*pi-pi/2;
% in(3)=in(1);
% in(3).shot = 160615026;
% in(3).line=2;
% in(3).color = {'g';'g'};%[66, 188, 244]./255};
% in(3).legend = [num2str(in(2).shot) ': -47.1kA (+)'];
% in(3).phaseShift = -pi/2;

% in(2).shot = 160615023;150625998;
% in(2).line =2; % line # NB: 1 is C III, 2 is O II, 3 is C III !
% in(2).legend = [num2str(in(1).shot) ' ' lines{in(1).line}];
% in(2).color = {'b';[66, 188, 244]./255};
% in(2).style = '-';
% in(2).error = 0; % 1 / 0 for errorbars
% in(2).velShift = -5; % SHIFT VELOCITY
% in(2).intScale = 1; % scale factor for intensity
% in(2).timeShift = 0; % ms, shift time base
% in(2).timeScale = 1e-3; % scale timebase to put into ms
% in(2).injTimeScale = 1;1e-3; % scale the injector time to ms
% in(2).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(2).doubleplot = [];[1:23; 24,26:47]; % plot coorespoinding impacts
% in(2).fftPlot = []; % FFT of signal, n frequencies
% in(3).shot = 160615024;150625998;
% in(3).line =2; % line # NB: 1 is C III, 2 is O II, 3 is C III !
% in(3).legend = [num2str(in(1).shot) ' ' lines{in(1).line}];
% in(3).color = {'g';[182, 244, 66]./255};
% in(3).style = '-';
% in(3).error = 0; % 1 / 0 for errorbars
% in(3).velShift = -5; % SHIFT VELOCITY
% in(3).intScale = 1; % scale factor for intensity
% in(3).timeShift = 0; % ms, shift time base
% in(3).timeScale = 1e-3; % scale timebase to put into ms
% in(3).injTimeScale = 1;1e-3; % scale the injector time to ms
% in(3).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(3).doubleplot = [];[1:23; 24,26:47]; % plot coorespoinding impacts
% in(3).fftPlot = []; % FFT of signal, n frequencies

% Highest performance HIT-SI3
% in(1).shot = 161018015;150625998;
% in(1).line =2; % line # NB: 1 is C III, 2 is O II, 3 is C III !
% in(1).legend = [num2str(in(1).shot) ' ' lines{in(1).line}];
% in(1).legend = [num2str(in(1).shot) ': +66kA'];
% in(1).color = {'r';'r'};%[225,105,0]./255};
% in(1).style = {'-','--'};
% in(1).error = 0; % 1 / 0 for errorbars
% in(1).velShift = -5; % SHIFT VELOCITY
% in(1).intScale = 2; % scale factor for intensity
% in(1).timeShift = 0; % ms, shift time base
% in(1).timeScale = 1e-3; % scale timebase to put into ms
% in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
% in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(1).doubleplot = [1];[1:23; 24,26:47]; % plot coorespoinding impacts
% in(1).fftPlot = [1]; % FFT of signal, n frequencies
% in(1).AnalysisTitle=[ 'HIT-SI3: 0-120-240 Phasing, Gain 2.9, ' lines{in(1).line}];
% in(1).phaseShift = 2*pi -pi/2;
% % 160728
% in(2)=in(1);
% in(2).shot = 161018013;
% in(2).line=2;
% in(2).color = {'b';'b'};%[66, 188, 244]./255};
% in(2).legend = [num2str(in(2).shot) ': +74.9'];
% in(2).phaseShift = -pi/2+pi;
% in(3)=in(1);
% in(3).shot = 160728013;
% in(3).color = {'g';'g'};%[182, 244, 66]./255};
% in(3).legend = [num2str(in(3).shot) ': +65.6kA'];


% % 129499 
% % note: need to change .*1e-6 to .*1e-3 in sinefit
in(1).shot = 129499;%150625998;
in(1).line = 3; % line # NB: 1 is C III, 2 is O II, 3 is C III !
in(1).legend = [num2str(in(1).shot) ': +90kA'];
in(1).color = {'r';'r'};%[225,105,0]./255};
in(1).style = {'-','--'};
in(1).error = 0; % 1 / 0 for errorbars
in(1).velShift = 0; % SHIFT VELOCITY
in(1).intScale = 1; % scale factor for intensity
in(1).timeShift = 0; % ms, shift time base
in(1).timeScale = 1;1e-3; % scale timebase to put into ms
in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
in(1).doubleplot = [];[1:23; 24,26:47]; % plot coorespoinding impacts
in(1).fftPlot = []; % FFT of signal, n frequencies
in(1).AnalysisTitle=['HIT-SI: 0-90 Phasing, ' lines{in(1).line+1}]; 
in(1).phaseShift =-pi/2;
in(1).error=0;
% 129450, 129451
% in(2) = in(1);
% in(2).shot=129450;
% in(2).line=1;
% in(2).legend = [num2str(in(2).shot) ': +79kA'];
% in(2).color = {'b';'b'};%[66, 188, 244]./255};
in(2)=in(1);
in(2).line=2
in(2).shot = 129496;
in(2).color = {'b';'b'};%[182, 244, 66]./255};
in(2).legend = [num2str(in(2).shot) ': -76kA'];
% Supression of Dissenting data:
phaseSupress(1:2,:)=22;
flowSupress(2,:)=22;
dispSupress(1:2,:)=22;
tempSupress(1:2,:)=[21,22;21,22];


%% Low Performance HIT-SI / Comparison
% in(1).shot = 129499;%150625998;
% in(1).line = 3; % line # NB: 1 is C III, 2 is O II, 3 is C III !
% in(1).legend = [num2str(in(1).shot) ': +90kA'];
% in(1).color = {'r';'r'};%[225,105,0]./255};
% in(1).style = {'-','--'};
% in(1).error = 0; % 1 / 0 for errorbars
% in(1).velShift = 5; % SHIFT VELOCITY
% in(1).intScale = 1; % scale factor for intensity
% in(1).timeShift = 0; % ms, shift time base
% in(1).timeScale = 1;1e-3; % scale timebase to put into ms
% in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
% in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(1).doubleplot = [];[1:23; 24,26:47]; % plot coorespoinding impacts
% in(1).fftPlot = [1]; % FFT of signal, n frequencies
% in(1).AnalysisTitle='HIT-SI: 0-90 Phasing, OII';
% in(1).phaseShift =-pi/2;
% % 129450, 129451
% in(2) = in(1);
% in(2).shot=129496;
% in(2).line=2;
% in(2).legend = [num2str(in(2).shot) ': -76kA'];
% in(2).color = {'b';'b'};%[66, 188, 244]./255};
% in(3)=in(1);
% in(3).shot = 129450;
% in(3).color = {'g';'g'};%[182, 244, 66]./255};
% in(3).legend = [num2str(in(3).shot) ': +78kA'];

% Low Performance HIT-SI3. 0-120-240
% in(1).shot = 151217024;150625998;
% in(1).line =2; % line # NB: 1 is C III, 2 is O II, 3 is C III !
% in(1).legend = [num2str(in(1).shot) ' ' lines{in(1).line}];
% in(1).legend = [num2str(in(1).shot) ': -41.1kA'];
% in(1).color = {'r';'r'};%[225,105,0]./255};
% in(1).style = {'-','--'};
% in(1).error = 0; % 1 / 0 for errorbars
% in(1).velShift = -5; % SHIFT VELOCITY
% in(1).intScale = 2; % scale factor for intensity
% in(1).timeShift = 0; % ms, shift time base
% in(1).timeScale = 1e-3; % scale timebase to put into ms
% in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
% in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(1).doubleplot = [1];[1:23; 24,26:47]; % plot coorespoinding impacts
% in(1).fftPlot = [1]; % FFT of signal, n frequencies
% in(1).AnalysisTitle=[ 'HIT-SI3: 0-120-240 Phasing, Gain 1.5, ' lines{in(1).line}];
% in(1).phaseShift = 2*pi -pi/2 + pi; % Add pi to make compairing neg-pos easier
% % 160728
% in(2)=in(1);
% in(2).shot = 151217025;
% in(2).line=2;
% in(2).color = {'b';'b'};%[66, 188, 244]./255};
% in(2).legend = [num2str(in(2).shot) ': -41.2kA'];
% in(2).phaseShift = -pi/2+pi;
% in(3)=in(1);
% in(3).shot = 151217026;
% in(3).phaseShift = 4*pi -pi/2;
% in(3).color = {'g';'g'};%[182, 244, 66]./255};
% in(3).legend = [num2str(in(3).shot) ': +38.1kA'];



% % in(2).shot = 151217026;%150625998;
% % in(2).line = 2; % line # NB: 1 is O II, 2 is C III !
% % in(2).legend = [num2str(in(1).shot) ' O II'];
% % in(2).color = {'b';[0,109,200]./255};
% % in(2).style = '-';
% % in(2).error = 0; % 1 / 0 for errorbars
% % in(2).velShift = 5; % SHIFT VELOCITY
% % in(2).intScale = 1; % scale factor for intensity
% % in(2).timeShift = 0; % ms, shift time base
% % in(2).timeScale = 1e-3; % scale timebase to put into ms
% % in(2).injTimeScale = 1e0;1e-3; % scale the injector time to ms
% % in(2).injScale = 1e0; 1e-3; % scale the inj current into kA
% % in(2).doubleplot = [1:23; 24,26:47]; % plot coorespoinding impacts
% % in(2).fftPlot = [1]; % FFT of signal, n frequencies

% in(3).shot = 151217026;%150625998;
% in(3).line = 3; % line # NB: 1 is O II, 2 is C III !
% in(3).legend = [num2str(in(1).shot) ' O II'];
% in(3).color = {[0,128,102]./255;[117,186,102]./255};
% in(3).style = '-';
% in(3).error = 0; % 1 / 0 for errorbars
% in(3).velShift = 5; % SHIFT VELOCITY
% in(3).intScale = 1; % scale factor for intensity
% in(3).timeShift = 0; % ms, shift time base
% in(3).timeScale = 1e-3; % scale timebase to put into ms
% in(3).injTimeScale = 1e0;1e-3; % scale the injector time to ms
% in(3).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(3).doubleplot = [1:23; 24,26:47]; % plot coorespoinding impacts
% in(3).fftPlot = [1]; % FFT of signal, n frequencies
% 
% in(2).shot = 12949610;
% in(2).line = 2; % line #
% in(2).legend = '129496 C III';
% in(2).color = 'r';
% in(2).style = '-';
% in(2).error = 1; % 1 / 0 for errorbars
% in(2).velShift = 0; % SHIFT VELOCITY
% in(2).intScale = 1; % scale factor for intensity
% in(2).timeShift = 0; % ms, shift time base
% in(1).shot = 160728011;150625998;
% in(1).line =1; % line # NB: 1 is C III, 2 is O II, 3 is C III !
% in(1).legend = [num2str(in(1).shot) ' ' lines{in(1).line}];
% in(1).legend = [num2str(in(1).shot) ': +65.5kA'];
% in(1).color = {'r';'r'};%[225,105,0]./255};
% in(1).style = {'-','--'};
% in(1).error = 0; % 1 / 0 for errorbars
% in(1).velShift = -5; % SHIFT VELOCITY
% in(1).intScale = 2; % scale factor for intensity
% in(1).timeShift = 0; % ms, shift time base
% in(1).timeScale = 1e-3; % scale timebase to put into ms
% in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
% in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(1).doubleplot = [1];[1:23; 24,26:47]; % plot coorespoinding impacts
% in(1).fftPlot = [1]; % FFT of signal, n frequencies
% in(1).AnalysisTitle=[ 'HIT-SI3: 0-120-240 Phasing, ' lines{in(1).line}];
% in(1).phaseShift = 2*pi -pi/2;
% % % 160728
% in(2)=in(1);
% in(2).shot = 160728012;
% in(2).line=1;
% in(2).color = {'b';'b'};%[66, 188, 244]./255};
% in(2).legend = [num2str(in(2).shot) ': -52.3kA'];
% in(2).phaseShift = -pi/2;
% in(3)=in(1);
% in(3).shot = 160728013;
% in(3).color = {'g';'g'};%[182, 244, 66]./255};
% in(3).legend = [num2str(in(3).shot) ': +65.6kA'];

%% 160609010 NIMROD
% in(1).shot = 8160609011;150625998;
% in(1).line =1; % line # NB: 1 is C III, 2 is O II, 3 is C III !
% in(1).legend = [num2str(in(1).shot) ' ' lines{in(1).line}];
% in(1).legend = [num2str(in(1).shot) ': +67.4kA'];
% in(1).color = {'r';'r'};%[225,105,0]./255};
% in(1).style = {'-','--'};
% in(1).error = 0; % 1 / 0 for errorbars
% in(1).velShift = -5; % SHIFT VELOCITY
% in(1).intScale = 2; % scale factor for intensity
% in(1).timeShift = 0; % ms, shift time base
% in(1).timeScale = 1;1e-3; % scale timebase to put into ms
% in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
% in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(1).doubleplot = [1];[1:23; 24,26:47]; % plot coorespoinding impacts
% in(1).fftPlot = [1]; % FFT of signal, n frequencies
% in(1).AnalysisTitle=[ 'NIMROD: 0-120-240 Phasing, D, Iso-Visc.'];
% in(1).phaseShift = 2*pi-pi/2;
% % 160728
% in(2)=in(1);
% in(2).shot = 8160609009;
% in(2).line=1;
% in(2).color = {'b';'b'};%[66, 188, 244]./255};
% in(2).legend = [num2str(in(2).shot) ': -52.3kA'];
% in(2).phaseShift = -pi/2;

%% 129499 NIMROD 
% in(1).shot = 8129499;%150625998;
% in(1).line = 1; % line # NB: 1 is C III, 2 is O II, 3 is C III !
% in(1).legend = [num2str(in(1).shot) ': +90kA'];
% in(1).color = {'r';[225,105,0]./255};
% in(1).style = {'-','--'};
% in(1).error = 0; % 1 / 0 for errorbars
% in(1).velShift = 5; % SHIFT VELOCITY
% in(1).intScale = 1; % scale factor for intensity
% in(1).timeShift = 0; % ms, shift time base
% in(1).timeScale = 1;1e-3; % scale timebase to put into ms
% in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
% in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
% in(1).doubleplot = [1];[1:23; 24,26:47]; % plot coorespoinding impacts
% in(1).fftPlot = [1]; % FFT of signal, n frequencies
% in(1).AnalysisTitle='HIT-SI: 0-90 Phasing, OII';
% in(1).phaseShift =2*pi-pi/2;
% 129450, 129451
% in(2) = in(1);
% in(2).shot=129450;
% in(2).line=1;
% in(2).legend = [num2str(in(2).shot) ': +89kA'];
% in(2).color = {'b';'b'};%[66, 188, 244]./255};
% in(3)=in(1);
% in(3).shot = 129496;
% in(3).color = {'g';'g'};%[182, 244, 66]./255};
% in(3).legend = [num2str(in(3).shot) ': -76kA'];



 plotType = 1; % Velocity
% plotType = 2; % Temperature
% plotType = 3; % Intensity
%plotType = 4; % Displacement

timebound = [1.1, 2.0]; % [ms]
% timebound = [1.2 2.3]; % [ms]

CutPow  = .1; % FFT Power Cuttoff Percentage. 


saving = 0;
plotFFT = 1;
plotCurrents = 1;
plotSanityPhase = 1;
plotAverages = 0;
compactCurrents = 1;
plotTor = 0; % The line plot will show the differenve between fibers
plotExplainReconst=0;% will plot the data and reconstruction from line n to a seperate plot
driveDirection=1*(plotType==1); % upper fiber dropped by -pi, phase now references "positive toroidal flow magnitude"
flipLoImpact=[47].*(in(1).shot>229499); % correct for lower fiber being rotated 180 about impact 30 
% the bracketed number is the pre-trimmed index of this impact. The
% following data is hardcoded for 160728017 data range
plotError = 1;
includeTemp = 1; % include temperature in the velocity plot

if isempty(in(1).fftPlot)
    Analysis=1;
else
    Analysis = 2; % Analyze torroidal flow, Amplitude, phasing. Replaces "Averages"
end


saveFile = ['T:\IDS\Analysis Repository\' num2str(in(1).shot)];

% Set these: ( add :2: if you only want every other line
% chan_range = 50:58;
% defaults
chan_ranget = [7:26];
chan_rangep = [43:62];
xlim = [0,50];
% 151217026 values:
if in(1).shot == 151217026
    chan_ranget = [13:26];[3:62];[8:26]; % toroidal, mohawk port in midplane
    chan_rangep = [49:62];[44:62]; % poloidal
    timebound=[.8,2.0];
elseif in(1).shot == 151217024
    chan_ranget = [10:26];
    chan_rangep = [46:62];
    timebound=[.8,2.0];
elseif in(1).shot == 151217025
    chan_ranget = [13:26];
    chan_rangep = [49:62];
    timebound=[.8,2.0];
elseif in(1).shot == 151217020 || in(1).shot == 151217021 || in(1).shot == 151217019 || in(1).shot == 151217022|| in(1).shot == 151217016
    chan_ranget = [15:25];
    chan_rangep = [51:61];
elseif in(1).shot <= 139499
    chan_ranget = [4:27];
    chan_rangep = []; % we dont want the axial port
    timebound = [1.6 2.0];
elseif in(1).shot ==8129499
    timebound = [0.3,0.66];
    chan_ranget = [1:36];
    chan_rangep = [37:72];
    xlim=[-20,60]
elseif in(1).shot == 160525016
    timebound = [1.5,2.0];
elseif in(1).shot == 160525017
    chan_ranget = [10:26];
    chan_rangep = [46:62];
    timebound = [1.35,2.0];
elseif in(1).shot >= 160526030 && in(1).shot <= 160526038
    timebound = [1.35,2.0];
    chan_ranget = [12:23];
    chan_rangep = [48:59];
elseif in(1).shot == 160615026
    timebound = [1.3,2.25];
    chan_ranget = [9:23];
    chan_rangep = [45:59];
elseif in(1).shot == 160615025
    timebound = [1.2,2.05];
    chan_ranget = [10:24];
    chan_rangep = [46:60];
elseif in(1).shot >= 160615022 && in(1).shot <= 160615024
    timebound = [1.75,2.05];
    chan_ranget = [10:24];
    chan_rangep = [46:60];
elseif in(1).shot == 160622021
    timebound = [1.2,2.05];
    chan_ranget = [10:24];
    chan_rangep = [46:60];
elseif in(1).shot >= 160728009 && in(1).shot <= 160728024
    timebound = [1.65,2.00];
%     chan_ranget = [8:24];
%     chan_rangep = [44:60];
    chan_ranget = [10:32];
    chan_rangep = [37:59];
elseif in(1).shot == 8160609009 %Nimrod
    timebound = [1.0,1.25];
    chan_ranget = [1:36];
    chan_rangep = [37:72];
    xlim=[-20,60]
elseif in(1).shot == 8160609010 %Nimrod
    timebound = [1.0,1.25];
    chan_ranget = [1:36];
    chan_rangep = [37:72];
    xlim=[-20,60]
elseif in(1).shot == 8160609011 %Nimrod
    timebound = [0.5,0.77];
    chan_ranget = [1:36];
    chan_rangep = [37:72];
    xlim=[-20,60]
elseif in(1).shot >= 161018010 && in(1).shot <= 161018023
    timebound = [1.6,2.0];
    chan_ranget = [7:17,19:29];
    chan_rangep = [43:65];
    xlim=[-20,60];
end
chan_range = [chan_ranget, chan_rangep];
%chan_range = chan_ranget;
%chan_range = chan_rangep
%chan_range = [1:68];
%chan_range = [8:2:24];
% chan_range = [43:53];

velSpace = 15; % km/s
intSpace = 20; % arb.u.
tempSpace = 20; % eV

%% set up figure
fntsz = 14;
lnwdth = 1.5;
errWdth = 500; % errorbar width setting


S = get(0, 'ScreenSize');
analysisHeight = S(4) - 1;
figureWidth = (S(3) - 12)/2.25;
Widen = 150; % 100 or zero usually
colors = [ 0 0 1; 12/255 117/255 0; 1 0 0 ];
h = figure('Visible', 'on', 'Name', ['MULTIPLOT-Lines: ' num2str(in(1).line)], 'Position',...
    [5, 35, figureWidth, analysisHeight], 'Color', [1 1 1]);
colorOrder = get(groot,'defaultAxesColorOrder');
partialColor = {'r',colorOrder(2,:); 'b',colorOrder(1,:);,'g',colorOrder(5,:)};
if ~compactCurrents
    ax = axes('Parent', h, 'Position', [0.075, 0.08, 0.8, 0.85]);
else
    ax = axes('Parent', h, 'Position', [0.075, 0.28, 0.8, 0.65]);
end

hold on;

if plotAverages || Analysis ==1
    h2 = figure('Visible', 'on', 'Name', ['MULTIPLOT-Toroidal Flow: ' num2str(in(1).line)], 'Position',...
        [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
    ax2 = axes('Parent', h2, 'Position', [0.075, 0.15, 0.85, 0.75]);
    hold on;
    grid on;
end
if Analysis==1 && ~isempty(in(1).fftPlot)
     h3 = figure('Visible', 'on', 'Name', ['MULTIPLOT-Phase: ' num2str(in(1).line)], 'Position',...
        [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
    ax3 = axes('Parent', h3, 'Position', [0.075, 0.15, 0.85, 0.75]); hold on; grid on;
    h4 = figure('Visible', 'on', 'Name', ['MULTIPLOT-Displacement: ' num2str(in(1).line)], 'Position',...
        [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
    ax4 = axes('Parent', h4, 'Position', [0.075, 0.15, 0.85, 0.75]); hold on; grid on;
     h5 = figure('Visible', 'on', 'Name', ['MULTIPLOT-FlowRanges: ' num2str(in(1).line)], 'Position',...
        [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
    ax5 = axes('Parent', h5, 'Position', [0.075, 0.15, 0.85, 0.75]);
    hold on;
    grid on;
elseif Analysis==2 && ~isempty(in(1).fftPlot)
    h2=figure;
    % Main Analysis Figure
    h6 = figure('Visible', 'on', 'Name', ['MULTIPLOT-Analysis: ' in(1).legend], 'Position',...
    [5, 1+25, figureWidth-450+Widen, analysisHeight-25], 'Color', [1 1 1]);
    ax6 = axes('Parent',h6,'Position',[0.15,.08,.8,.25]); hold on; grid on;box on;
    ax7 = axes('Parent',h6,'Position',[0.15,.38,.8,.25]); hold on; grid on;box on;
    ax8 = axes('Parent',h6,'Position',[0.15,.68,.8,.25]); hold on; grid on;box on;
    %txt=text(10,550,['Analysis: ' in(1).legend],'fontweight','bold','fontsize',13);
    mBox=uicontrol('Style','Text');% Make the title
    set(mBox,'String',in(1).AnalysisTitle); set(mBox,'fontweight','bold'); set(mBox,'fontsize',13);
    set(mBox,'Position', [ 30+ceil(Widen/2),960,400,20]);set(mBox,'BackgroundColor',[1 1 1]);
    if plotType==1
        title(ax8,'Toroidal Displacement Phase');
        title(ax7,'Toroidal Flow');
        title(ax6,'Maximum Displacement');
    elseif plotType==2
        title(ax8,'Temperature Phase');
        title(ax7,'Temperature Profile');
        title(ax6,'Temperature Oscillations');
    end
    if includeTemp
    ax17=axes('Parent',h6,'Position',[0.15,.08,.8,.18]); hold on; grid on;box on;
    title(ax17,'Temperature Profile');
    set(ax6,'Position',[0.15,.30,.8,.18]); 
    set(ax7,'Position',[0.15,.52,.8,.18]);
    set(ax8,'Position',[0.15,.74,.8,.18]);
    end
    
    % Plot FFT Power Spectrum
    if plotFFT
    h7 = figure('Visible', 'on', 'Name', ['MULTIPLOT-FFT Spectrum: ' num2str(in(1).line)], 'Position',...
        [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
    ax9 = axes('Parent', h7, 'Position', [0.075, 0.15, 0.85, 0.75]); hold on; grid on;box on;
    end
    
    if plotSanityPhase
    % Show the phases as calculated various methods
    h10=figure('Visible', 'on', 'Name', ['Phase: Sanity Check ' in(1).legend], 'Position',...
    [5, 1, figureWidth-300, analysisHeight], 'Color', [1 1 1]);
    ax10 = axes('Parent',h10,'Position',[0.15,.08,.8,.25]); hold on; grid on;box on;
    ax11 = axes('Parent',h10,'Position',[0.15,.38,.8,.25]); hold on; grid on;box on;
    ax12 = axes('Parent',h10,'Position',[0.15,.68,.8,.25]); hold on; grid on;box on;
    % Reconstruct a few lines as a sanity check    
    h11=figure('Visible', 'on', 'Name', ['Phase: Reconstruction ' in(1).legend], 'Position',...
    [5, 1, figureWidth-300, analysisHeight], 'Color', [1 1 1]);
    ax13 = axes('Parent',h11,'Position',[0.15,.08,.8,.25]); hold on; grid on;box on;
    ax14 = axes('Parent',h11,'Position',[0.15,.38,.8,.25]); hold on; grid on;box on;
    ax15 = axes('Parent',h11,'Position',[0.15,.68,.8,.25]); hold on; grid on;box on;
    end
    
    if plotExplainReconst
    % Data Reconstruction Figure
    h12=figure('Visible','on','Name','Explainatory Reconstruction');
    ax16 = axes('Parent',h12); hold on; grid on; box on;
    end
%     h3=figure;
%     h4=figure;
%     h5=figure;
end
figure(h); % make first figure current

%% Load and Plot Data
clear param
for n = 1:length(in)
    n
    saveDat(n).title= in(n).AnalysisTitle;
    saveDat(n).shot=in(n).shot;
    clear data zeroline time
    %addpath('T:\IDS\Data Repository');
    
    if ~(exist('dat','var') && strcmp(dat(1).title,['Shot ' num2str(in(n).shot)])) % speed up runtime if data is already loaded
        input=load(['T:\IDS\Data Repository\dat' num2str(in(n).shot) '10.mat']); % Real HIT-SI Data
        dat=input.dat;
        if flipLoImpact ~=0 % flip the lower array about its centeral impact.
               breakInd = 30; % Assume that the fiber break occurs at index thirty
              impacts=importdata('impacts5.mat')';
                dat(in(n).line).vel(:,breakInd:end) = [dat(in(n).line).vel(:,end:-1:end-2*(length(dat(1).impacts)-flipLoImpact))...
                    dat(in(n).line).vel(:,length(dat(1).impacts)-flipLoImpact+breakInd:-1:breakInd) ];
                dat(in(n).line).temp(:,breakInd:end) = [dat(in(n).line).temp(:,end:-1:end-2*(length(dat(1).impacts)-flipLoImpact)) ...
                    dat(in(n).line).temp(:,length(dat(1).impacts)-flipLoImpact+breakInd:-1:breakInd) ];
                 dat(in(n).line).velU(:,breakInd:end) = [dat(in(n).line).velU(:,end:-1:end-2*(length(dat(1).impacts)-flipLoImpact))...
                    dat(in(n).line).velU(:,length(dat(1).impacts)-flipLoImpact+breakInd:-1:breakInd) ];
                dat(in(n).line).tempU(:,breakInd:end) = [dat(in(n).line).tempU(:,end:-1:end-2*(length(dat(1).impacts)-flipLoImpact)) ...
                    dat(in(n).line).tempU(:,length(dat(1).impacts)-flipLoImpact+breakInd:-1:breakInd) ];
                dat(in(n).line).velL(:,breakInd:end) = [dat(in(n).line).velL(:,end:-1:end-2*(length(dat(1).impacts)-flipLoImpact))...
                    dat(in(n).line).velL(:,length(dat(1).impacts)-flipLoImpact+breakInd:-1:breakInd) ];
                dat(in(n).line).tempL(:,breakInd:end) = [dat(in(n).line).tempL(:,end:-1:end-2*(length(dat(1).impacts)-flipLoImpact)) ...
                    dat(in(n).line).tempL(:,length(dat(1).impacts)-flipLoImpact+breakInd:-1:breakInd) ];
                dat(1).impacts(breakInd:end) = [dat(1).impacts(end-2*(length(dat(1).impacts)-flipLoImpact):end);...
                    impacts(63:62+(-length(dat(1).impacts)+2*flipLoImpact-breakInd))];
         end

        dat = trimRange(dat, chan_range, plotError,timebound.*(1./in(n).timeScale),[]); % for some reason, this wont save to workspace
        assignin('base','dat',dat);
        %try;saveDat(n).VelError = dat(in(n).line).velU;end
        %try;saveDat(n).TempError = dat(in(n).line).tempU;end
    end
    temp = dat(1).vel;
    
    Itor(:,n) = dat(1).Itor;
    
    %switch plotType
        if plotType==1 || plotType==4 || ( plotType ==2 && ~isempty(in(1).fftPlot))
            if ~isempty(in(n).doubleplot) && isempty(in(1).fftPlot)
                %data(1:length(dat(1).time),:) = dat(in(n).line).vel(:,in(n).doubleplot(1,:));
                %data(length(dat(1).time)+1:2*length(dat(1).time),:) = ...
                %    dat(in(n).line).vel(:,in(n).doubleplot(2,:));
                if plotType ==1
                data(1:length(dat(1).time),:) = dat(in(n).line).vel(:,1:(length(dat(1).impacts))/2)+in(n).velShift;
                data(length(dat(1).time)+1:2*length(dat(1).time),:) = ...
                   dat(in(n).line).vel(:,(length(dat(1).impacts)/2)+1:end)+in(n).velShift;
                elseif plotType==4
                    data(1:length(dat(1).time),:) = dat(in(n).line).vel(:,1:(length(dat(1).impacts))/2);
                end
            elseif ~isempty(in(1).fftPlot)
                dat(in(n).line).vel = averageNans(dat(in(n).line).vel)+in(n).velShift; % remove nans
                dat(in(n).line).temp = averageNans(dat(in(n).line).temp); % remove nans
                saveDat(n).rawVel = dat(in(n).line).vel;
                saveDat(n).rawTemp = dat(in(n).line).temp;
               
                if ~isempty(in(n).doubleplot)
                    % Find where each fiber bundle begins and ends.
                    doubleplot(1,:) = 1:(length(dat(1).impacts))/2;
                    doubleplot(2,:) = (length(dat(1).impacts)/2)+1:length(dat(1).impacts);

                    % initialize Sine_Fit parameters
                    param(:,1,n) = dat(1).impacts(doubleplot(1,:));
                    param(:,6,n) = dat(1).impacts(doubleplot(2,:));
                else
                    doubleplot(1,:) = 1:length(dat(1).impacts);
                    param(:,1,n) = dat(1).impacts(doubleplot(1,:));
                end
%                 if flipLoImpact == 1 % flip the lower array about its centeral impact.
%                     dat(in(n).line).vel(:,doubleplot(2,:)) = dat(in(n).line).vel(:,doubleplot(2,end):-1:doubleplot(2,1));
%                     dat(in(n).line).temp(:,doubleplot(2,:)) = dat(in(n).line).temp(:,doubleplot(2,end):-1:doubleplot(2,1));
%                 end
                 % Param: impacts offset amplitude phase, frequency
                 if ~exist('param','var')
                        param = zeros(length(doubleplot),10,length(in));
                 end
                display('Computing FFT');
                pRel = zeros(length(doubleplot),2);
                dPar = zeros(length(in),size(doubleplot,2),size(doubleplot,2),4);
                for i = 1:length(doubleplot) % Loop through impacts ( fits both arrays)
                    % using fftf method
%                    offset = mean(dat(in(n).line).vel(:,in(n).doubleplot(1,i)));
%                    [data(1:length(dat(1).time),i), f, y, y2] = ...
%                        fftf(dat(1).time.*1e-6, dat(in(n).line).vel( ...
%                        :,in(n).doubleplot(1,i)) - offset, 50e3,3,10e3,1);
%                    data(1:length(dat(1).time),i) = ...
%                    data(1:length(dat(1).time),i).* abs(max(y2)); % make real
%                    data(1:length(dat(1).time),i) = ...
%                        data(1:length(dat(1).time),i) + offset; % replace offset
%                    pause(1);
%                    offset = mean(dat(in(n).line).vel(:,in(n).doubleplot(2,i)));
%                    [data(length(dat(1).time)+1:2*length(dat(1).time),i), f, y, y2] = ...
%                        fftf(dat(1).time.*1e-6, dat(in(n).line).vel( ...
%                        :,in(n).doubleplot(2,i))- offset, 50e3,3,10e3,1);
%                    data(length(dat(1).time)+1:2*length(dat(1).time),i) = ...
%                        data(length(dat(1).time)+1:2*length(dat(1).time),i).* abs(max(y2));
%                    data(length(dat(1).time)+1:2*length(dat(1).time),i) = ...
%                        data(length(dat(1).time)+1:2*length(dat(1).time),i) + offset;
%                    pause(1);
                    % using SineFit method
                    %try
                    
                    %%%%%% Upper Array %%%%%%%%%%%
                    % Calculate sine fit for uppe array
                    size(dat(in(n).line).vel);
                        if plotType==1
                            signal = dat(in(n).line).vel(:,doubleplot(1,i));
                        elseif plotType == 2
                            signal = dat(in(n).line).temp(:,doubleplot(1,i));
                        end
                        if exist('savSig','var');size(savSig)
                        end
                        size(signal)
                        %savSig(1,i,n,:)=signal;
                        Fsamp = 1/(mean(diff(dat(1).time.*(in(n).timeScale.*1e-3))));
                        offset = nanmean(signal);
                        amp = max(signal)-offset;
                        freq = 14500;
                        phase = pi/2;
                        
                        
                        xfft(:,i,n,1)=fft(signal-offset);
                        P1=xfft(1:floor(length(signal)/2)+1,i,n,1);
                        P2=abs(P1);
                        P2(2:end-1)=2*P2(2:end-1);
                        f=Fsamp*(0:length(signal)-1)/length(signal);
                        [Y,I]=max(P2);
                        amp = P2(I)/length(signal);
                        phase = mod(unwrap(angle(P1(I))),2*pi);
                        data_fft(1:length(dat(1).time),i,n)=offset+amp*sin(f(I)*2*pi*dat(1).time.*(in(n).timeScale.*1e-3) + phase);
                        param_fft(i,2:5,n) = [offset, amp, phase, f(I)];
                        if f(I)>20000 ; f(I)=f(I)/2; end % hit harmoinc
                        
                        guess(i,:,n,1) = [offset,max(signal)-offset,pi/2,freq];
                        
                        [param(i,2:5,n),data(1:length(dat(1).time),i)] = sine_fit( ...
                            dat(1).time'.*(in(n).timeScale.*1e-3),signal',[nan,nan,nan,freq], ...
                            guess(i,:,n,1),0);
                        % Use FFT Result
                        
                        if param(i,3,n)<0 % 180degree phase
                            param(i,3,n)=-param(i,3,n);
                            param(i,4,n)=param(i,4,n)+pi;
                            disp([' WARNING: NEGATIVE AMPLITUDE @ n=' num2str(n) ', Line 1, Impact: ' num2str(dat(1).impacts(i))]);
                        end
                        param(i,3,n)= param_fft(i,3,n);
                        
                        
                        % Calculate fft reconstruction validity
                        harm1 = bandpower(signal-offset,Fsamp,[13,17].*1e3);
                        harm2 = bandpower(signal-offset,Fsamp,[26,31].*1e3);
                        harm3=0;harm4=0;
                        try harm3 = bandpower(signal-offset,Fsamp,[41,46].*1e3);end
                        try harm4 = bandpower(signal-offset,Fsamp,[55.5,60.5].*1e3);end
                        %harm5 = bandpower(signal-offset,Fsamp,[70,75].*1e3);
                        ptot =  bandpower(signal-offset,Fsamp,[0,Fsamp*(length(signal)-2)/(2*length(signal)) ]); % Nyquist
                        pRel(i,1) = (harm1+harm2+harm3+harm4)/ptot;
                        % calculate sine fit from FFT (again)
                        
                        % if the fit isnt valid, dont plot it
                        try data(1:length(dat(1).time),i.*(pRel(i,1)<CutPow))=signal;end
                        
                        % attempt lm error analysis
%                         dp = [0.001, 0.001, 0.001, 0.001]; % fractional increment of 'p' for numerical derivatives
%                         [p_fit, Chi_sq, dPar(n,i,1,:), ~, corr, R2, cvg_hst] = ...
%                         lm(@SineFitLM, param(i,2:5,n), dat(1).time.*(in(n).timeScale.*1e-3), signal', 0.001, dp);%, p_min,p_max,0)
                        %pause(1);
                        
                        % Calculate RMS Error
                        % Calculate RMS Errorfor rmsPhase = 1:200
                             for rmsPhase = 1:200
                                 %RMS(2,i,n,rmsPhase) = sqrt(mean( (signal'-(param(i,2,n)+param(i,3,n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,4,n) + (-pi +pi*rmsPhase/10)))).^2));
                                 RMS(1,i,n,rmsPhase) = sqrt(mean( (signal'-(param(i,2,n)+param(i,3,n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,4,n) + (-pi +pi*rmsPhase/100))) ).^2 ));
                                 RMS_ideal(1,i,n,rmsPhase)= sqrt(mean( ((param(i,2,n)+param(i,3,n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,4,n)))-(param(i,2,n)+param(i,3,n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,4,n) + (-pi +pi*rmsPhase/100)) ) ).^2 ));
                             end
%                             [p_fit_dat(n,i,2,:), Chi_sq, dPar_dat(n,i,2,:), ~, corr, R2, cvg_hst] = ...
%                             lm(@singletGauss1DLM, [-(RMS(2,i,n,1)-RMS(2,i,n,10))*pi,0,1.5,RMS(2,i,n,1)], (-pi +pi*(1:20)/10)', RMS(2,i,n,:), 0.0001, ones(1,4).*.001);%, p_min,p_max,0)
%                              p_fit_dat(n,i,1,:)=nlinfit((-pi +pi*(1:200)/100),squeeze(RMS(1,i,n,:))',@singletGauss1D,[-(RMS(1,i,n,1)-RMS(1,i,n,50))*pi,.5,1.5,RMS(1,i,n,1)]);
%                              [Y,I] = min( (squeeze(RMS_ideal(1,i,n,:))- RMS_ideal(1,i,n,1)/2).^2);% Find HWHM
%                              SigDev(n,i,1) = p_fit_dat(n,i,1,3)-abs(-pi +pi*(I)/100);
                             [p_fit_dat(n,i,1,:),R(n,i,1,:),J,COVB(n,i,1,:,:),MSE(n,i,1),ERRORMODELINFO(n,i,1)]=nlinfit((-pi +pi*(1:200)/100),squeeze(RMS(1,i,n,:))',@singletGauss1D,[-(RMS(1,i,n,1)-RMS(1,i,n,100))*pi,.5,1.5,RMS(1,i,n,1)]);
                             [p_fit_dat_ideal(n,i,1,:),R_ideal(n,i,1,:),J,COVB_ideal(n,i,1,:,:),MSE_ideal(n,i,1),ERRORMODELINFO_ideal(n,i,1)]=nlinfit((-pi +pi*(1:200)/100),squeeze(RMS_ideal(1,i,n,:))',@singletGauss1D,[-(RMS_ideal(1,i,n,1)-RMS_ideal(1,i,n,100))*pi,.5,1.5,RMS_ideal(1,i,n,1)]);
                             
%                              [Y,I] = min( (squeeze(RMS_ideal(2,i,n,:))- RMS_ideal(2,i,n,1)/2).^2);% Find HWHM
%                              SigDev(n,i,2) = p_fit_dat(n,i,2,3)-abs(-pi +pi*(I)/100);
                                SigDev(n,i,1) = abs(p_fit_dat(n,i,1,3)-p_fit_dat_ideal(n,i,1,3)); % Delta sigmas is error
                        %%%%%% Lower Array %%%%%%%%%%%%%%%
                        if in(n).doubleplot
                            % Fit Sine
                            if plotType==1
                                signal = dat(in(n).line).vel(:,doubleplot(2,i));
                            elseif plotType == 2
                                signal = dat(in(n).line).temp(:,doubleplot(2,i));
                            end
                            %savSig(2,i,n,:)=signal;
                            offset = mean(signal);
                            amp = max(signal)-offset;
                            
                            xfft(:,i,n,2)=fft(signal-offset);
                            P1=xfft(1:floor(length(signal)/2)+1,i,n,2);
                            P2=abs(P1);
                            P2(2:end-1)=2*P2(2:end-1);
                            f=Fsamp*(0:length(signal)-1)/length(signal);
                            [Y,I]=max(P2);
                            amp = P2(I)/length(signal);
                            phase = mod(unwrap(angle(P1(I))),2*pi);
                            data_fft(1:length(dat(1).time),i,n)=offset+amp*sin(f(I)*2*pi*dat(1).time.*(in(n).timeScale.*1e-3) + phase);
                            param_fft(i,7:10,n) = [offset, amp, phase, f(I)];
                            if f(I)>20000 ; f(I)=f(I)/2; end % hit harmoinc
                            
                            guess(i,:,n,2) = [offset,max(signal)-offset,pi/2,freq];
                            [param(i,7:10,n),data(length(dat(1).time)+1:2*length(dat(1).time),i)] = ...
                                sine_fit(dat(1).time'.*(in(n).timeScale.*1e-3),signal',[nan,nan,nan,freq], ...
                                guess(i,:,n,2),0);
                            %pause(1);
                            % Use the FFT Results
                            %param(i,8,n)= abs(param_fft(i,8,n));
                            if param(i,8,n)<0 % 180degree phase
                                param(i,8,n)=-param(i,8,n);
                                param(i,9,n)=param(i,9,n)+pi;
                                disp([' WARNING: NEGATIVE AMPLITUDE @ n=' num2str(n) ', Line 2, Impact: ' num2str(dat(1).impacts(i))]);
                            end
                            param(i,8,n)= (param_fft(i,8,n)); % needs to go here, otherwise sine of amp cant get fixed
                            % Calculate fft reconstruction validity
                            harm1 = bandpower(signal-offset,Fsamp,[13,17].*1e3);
                            harm2 = bandpower(signal-offset,Fsamp,[26,31].*1e3);
                            harm3=0;harm4=0;
                            try harm3 = bandpower(signal-offset,Fsamp,[41,46].*1e3);end
                            try harm4 = bandpower(signal-offset,Fsamp,[55.5,60.5].*1e3);end
                            ptot =  bandpower(signal-offset,Fsamp,[0,Fsamp*(length(signal)-2)/(2*length(signal))]);
                            pRel(i,2) = (harm1+harm2+harm3+harm4)/ptot;
                            xfft(:,i,n,2)=abs(fft(signal-offset))/length(signal);
                            % if the fit isnt valid, dont plot it
                            try data(length(dat(1).time)+1:2*length(dat(1).time)...
                                    ,i.*(pRel(i,2)<CutPow))=signal;end
                            
                            % attempt lm error analysis
%                             dp = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001]; % fractional increment of 'p' for numerical derivatives
%                             [p_fit, Chi_sq, dPar(n,i,1,:), ~, corr, R2, cvg_hst] = ...
%                             lm(@SineFitLM,  param(i,7:10,n), dat(1).time.*(in(n).timeScale.*1e-3), signal, 0.0001, dp);%, p_min,p_max,0)
%                              
                            % Calculate RMS Errorfor rmsPhase = 1:200
                             for rmsPhase = 1:200
                                 %RMS(2,i,n,rmsPhase) = sqrt(mean( (signal'-(param(i,2,n)+param(i,3,n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,4,n) + (-pi +pi*rmsPhase/10)))).^2));
                                 RMS(2,i,n,rmsPhase) = sqrt(mean( (signal'-(param(i,7,n)+param(i,8,n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,9,n) + (-pi +pi*rmsPhase/100))) ).^2 ));
                                 RMS_ideal(2,i,n,rmsPhase)= sqrt(mean( ((param(i,7,n)+param(i,8,n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,9,n)))-(param(i,7,n)+param(i,8,n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,9,n) + (-pi +pi*rmsPhase/100)) ) ).^2 ));
                             end
%                             [p_fit_dat(n,i,2,:), Chi_sq, dPar_dat(n,i,2,:), ~, corr, R2, cvg_hst] = ...
%                             lm(@singletGauss1DLM, [-(RMS(2,i,n,1)-RMS(2,i,n,10))*pi,0,1.5,RMS(2,i,n,1)], (-pi +pi*(1:20)/10)', RMS(2,i,n,:), 0.0001, ones(1,4).*.001);%, p_min,p_max,0)
                             [p_fit_dat(n,i,2,:),R(n,i,2,:),J,COVB(n,i,2,:,:),MSE(n,i,2),ERRORMODELINFO(n,i,2)]=nlinfit((-pi +pi*(1:200)/100),squeeze(RMS(2,i,n,:))',@singletGauss1D,[-(RMS(2,i,n,1)-RMS(2,i,n,100))*pi,.5,1.5,RMS(2,i,n,1)]);
                             [p_fit_dat_ideal(n,i,2,:),R_ideal(n,i,2,:),J,COVB_ideal(n,i,2,:,:),MSE_ideal(n,i,2),ERRORMODELINFO_ideal(n,i,2)]=nlinfit((-pi +pi*(1:200)/100),squeeze(RMS_ideal(2,i,n,:))',@singletGauss1D,[-(RMS_ideal(2,i,n,1)-RMS_ideal(2,i,n,100))*pi,.5,1.5,RMS_ideal(2,i,n,1)]);
                             
%                              [Y,I] = min( (squeeze(RMS_ideal(2,i,n,:))- RMS_ideal(2,i,n,1)/2).^2);% Find HWHM
%                              SigDev(n,i,2) = p_fit_dat(n,i,2,3)-abs(-pi +pi*(I)/100);
                                SigDev(n,i,2) = abs(p_fit_dat(n,i,2,3)-p_fit_dat_ideal(n,i,2,3)); % Delta sigmas is error
                        end
                            
                    
                    %catch
                    %    display(['halted at ' num2str(i)]);
                    %end
                end
                %%%%%%%% Find Phase of X/A Injector
                if n==1
                clear signal
                if (in(1).shot <129500)|| (in(1).shot == 8129499)
                    [Y,I1]=min( (dat(1).iinjxTime - timebound(1) ).^2 );
                    [Y,I2]=min( (dat(1).iinjxTime - timebound(2) ).^2 );
                    signal(:,1) = dat(1).iinjx(I1:I2);
                    signal(:,2) = dat(1).iinjy(I1:I2);
                    for i = 1:2
                        offset = mean(signal(:,i));
                        amp = max(signal(:,i))-offset;
                        [injParam(i,:),~] = ...
                            sine_fit(double(dat(1).iinjxTime(I1:I2))'.*1e-3,double(signal(:,i))',[nan,nan,nan,freq], ...
                            [offset,amp,0,freq],0);
                    end
                else
                    [Y,I1]=min( (dat(1).iinjaTime - timebound(1) ).^2 );
                    [Y,I2]=min( (dat(1).iinjaTime - timebound(2) ).^2 );
                    signal(:,1) = dat(1).iinja(I1:I2);
                    signal(:,2) = dat(1).iinjb(I1:I2);
                    signal(:,3) = dat(1).iinjc(I1:I2);
                    offset = mean(signal);
                    amp = max(signal)-offset;
                    for i = 1:3
                        offset = mean(signal(:,i));
                        amp = max(signal(:,i))-offset;
                        [injParam(i,:),~] = ...
                            sine_fit(double(dat(1).iinjaTime(I1:I2))'.*1e-3,double(signal(:,i))',[nan,nan,nan,freq], ...
                            [offset,amp,0,freq],0);
                    end
                end
                for i = 1:size(injParam,1)
                    if injParam(i,2)<0 % 180degree phase
                        injParam(i,2)=-injParam(i,2);
                        injParam(i,3)=injParam(i,3)+pi;
                        disp(['WARNING: INJECTOR ' num2str(i) ' NEGATIVE AMPLITUDE']);
                    end
                end
                
                end
                
            else % Just plot Lines
                if plotType==1
                    data = dat(in(n).line).vel + in(n).velShift;
                elseif plotType==4
                    for i = 1:size(dat(1).vel,2)
                    data(:,i) = cumtrapz(dat(1).time',dat(in(n).line).vel(:,i)-...
                        mean(dat(in(n).line).vel(:,i))).*100;
                    end
                end
            end
            titles = 'Velocities';
            offset = velSpace;
            units = 'km/s';
            sidebar = [num2str(offset) ' km/s per division'];
            if in(n).error
                errorL = dat(in(n).line).velL;
                errorU = dat(in(n).line).velU;
            end
        elseif plotType==2
            % Check Double Plotting
             if ~isempty(in(n).doubleplot) && isempty(in(1).fftPlot)
                doubleplot(1,:) = 1:(length(dat(1).impacts))/2;
                doubleplot(2,:) = (length(dat(1).impacts)/2)+1:length(dat(1).impacts);
                data(1:length(dat(1).time),:) = dat(in(n).line).temp(:,doubleplot(1,:));
                data(length(dat(1).time)+1:2*length(dat(1).time),:) = ...
                    dat(in(n).line).temp(:,doubleplot(2,:));
             else
            data = dat(in(n).line).temp;
             end
             
            titles = 'Temperatures';
            offset = tempSpace;
            units = 'eV';
            sidebar = [num2str(offset) ' eV per division'];
            if in(n).error
                errorL = dat(in(n).line).tempL;
                errorU = dat(in(n).line).tempU;
            end
        elseif plotType==3
            data = in(n).intScale * dat(in(n).line).int;
            titles = 'Intensities';
            sidebar = ['Arb.'];
            offset = intSpace;
            units = 'Arb.';
            if in(n).error
                errorL = dat(in(n).line).intL;
                errorU = dat(in(n).line).intU;
            end
        end
%     size(data)
    
    
    
    if in(n).doubleplot
         time(1:length(dat(1).time)) = dat(1).time.*in(1).timeScale + in(n).timeShift;
         time(length(dat(1).time)+1:2*length(dat(1).time)) =  time(1:length(dat(1).time));
%              dat(1).time.*in(1).timeScale + in(n).timeShift;
    else
        time = dat(1).time.*in(1).timeScale + in(n).timeShift;
    end

    %% Averages
    if plotAverages
        figure(h2) % make current
        
        % find index corresponding to time bounds
        nTimeLim(1) = find(dat(1).time >= timebound(1), 1);
        nTimeLim(2) = find(dat(1).time <= timebound(end), 1, 'last');
%         nTimeLim(2) = nTimeLim(2) - 1; % the above command goes one too far
        
        % Calculate all data
        for m = 1:size(data, 2);
            selection = data(nTimeLim(1):nTimeLim(2), m);
            dataAvg(m) = mean(selection(~isnan(selection)));
            dataStd(m) = std(selection(~isnan(selection)));
        end
        
        % Plot Data
        t2(n) = errorbar(dat(1).impacts, dataAvg, dataStd, 'color', in(n).color, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
        
        figure(h) % make other current
    end
    
    %% Analysis
    if Analysis
        %%%%%%%%%%%%%%%%%%%% Calculating %%%%%%%%%%%%%%%%%%%
        if plotType ==1 || ( plotType ==2 && ~isempty(in(1).fftPlot))
            if n==1
                figure(h2)
            end
            lines=1;
            if ~isempty(in(n).doubleplot)
                lines=2;
            end
            

            for i = 1:lines
                if in(n).fftPlot % if computing the fft, use offset, amp 
                    % Toroidal Current
                    dataAvg(:,i) = param(:,2+5*(i-1),n);
                    dataStd(:,i) = param(:,3+5*(i-1),n);
                    
                    %% Phase Modification Stuff
                    dataPhase(:,i) = param(:,4+5*(i-1),n);
 
                      % mod the phase
                    % Sanity check
                    if plotSanityPhase==1;eval(['plot(ax' num2str(9+n) ',dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,''color'',partialColor{n,i},''linewidth'',lnwdth)']);end
                    dataPhase(:,i) = mod(dataPhase(:,i),2*pi);
                    if plotSanityPhase==1;eval(['plot(ax' num2str(9+n) ',dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,''color'',partialColor{n,i},''linewidth'',lnwdth,''linestyle'','':'')']);end
                    % test automatic phase finding
                    for j = 1: size(param,1) % for all impacts
                        dataPhase(j,i) = find_Phase(param(j,1+(5*(i-1)): 5 + 5*(i-1),n));
                    end
                    if plotSanityPhase==1;eval(['plot(ax' num2str(9+n) ',dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,''color'',partialColor{n,i},''linewidth'',lnwdth,''linestyle'',''--'')']);end
                   % shift lines if a periodicity jump occurs
                    for j = 2:size(data, 2)
                        if (dataPhase(j,i) - dataPhase(j-1,i))>(pi)
                            dataPhase(j,i) = dataPhase(j,i)-(2*pi);
                            disp(['JUMP: -2Pi, Line: ' num2str(i) ' Impact: ' num2str(dat(1).impacts(j))]);
                        elseif (dataPhase(j,i) - dataPhase(j-1,i))<(-pi)
                           disp(['JUMP: +2Pi, Line: ' num2str(i) ' Impact: ' num2str(dat(1).impacts(j))]);
                            dataPhase(j,i) = dataPhase(j,i)+(2*pi);
                        end
                    end
                    if plotSanityPhase==1;eval(['plot(ax' num2str(9+n) ',dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,''color'',partialColor{n,i},''linewidth'',lnwdth,''linestyle'',''-.'')']);end
                    if plotSanityPhase==1;
                    legend(ax10,{'Raw Phase';'Mod 2\pi Phase';'Find\_Phase';'Jump Corrected';},'location','SouthWest');
                    set(ax11,'xticklabel',[]);
                    set(ax12,'xticklabel',[]);
                    eval(['ylabel(ax' num2str(9+n) ','' Phase [deg]'')']);
                    xlabel(ax10,'Impacts [cm]');
                   mBox1=uicontrol(h10,'Style','Text');% Make the title
                    set(mBox1,'String',['Phase Calculation Methods']); 
                    set(mBox1,'fontweight','bold'); set(mBox1,'fontsize',13);
                    set(mBox1,'Position', [ 80,960,450,20]);set(mBox1,'BackgroundColor',[1 1 1]);
                    eval(['title(ax' num2str(9+n) ', ''Shot: ' num2str(in(n).shot) ''')']);
                    end
                    dataPhase(:,i) = mod(param(:,4+5*(i-1),n),2*pi); % 2Pi Periodicity
                    
                    % Plot Reconstruction
                    if plotSanityPhase==1;
                    if in(1).shot >129500
                        %if i==2 % only plot the second time around
                            plotChan=9;
                        for j = 1:(n==plotExplainReconst)*1+1 % only plot the special explainatory reconstruction once
                            reconTime = linspace(dat(1).time(1).*in(n).timeScale,dat(1).time(end).*in(n).timeScale,length(dat(1).time)*10);
                            eval(['plot(ax' num2str(12+n+(j==2)*1) ',dat(1).time.*in(n).timeScale,dat(in(n).line).vel(:,plotChan),''color'',partialColor{n,1},''linewidth'',lnwdth)']);
                            eval(['plot(ax' num2str(12+n+(j==2)*1) ',dat(1).time.*in(n).timeScale,dat(in(n).line).vel(:,plotChan+(length(dat(1).impacts))/2),''color'',partialColor{n,2},''linewidth'',lnwdth,''linestyle'',''-'')']);
                            eval(['plot(ax' num2str(12+n+(j==2)*1) ',reconTime,param(plotChan,2+5*(i-1),n)+param(plotChan,3+5*(i-1),n)*sin(2*pi*reconTime.*1e-3.*14500 + dataPhase(plotChan,i)),''color'',partialColor{n,i},''linewidth'',lnwdth,''linestyle'',''--'')']);
                            %eval(['plot(ax' num2str(12+n) ',dat(1).time.*in(n).timeScale,param(plotChan,7,n)+param(plotChan,8,n)*sin(2*pi*dat(1).time.*1e-3.*in(n).timeScale.*14500 + dataPhase(plotChan,2)),''color'',partialColor{n,2},''linewidth'',lnwdth,''linestyle'',''--'')']);
                            eval(['title(ax' num2str(12+n) ', ''Shot: ' num2str(in(n).shot) ''')']);
                            set(ax15,'xticklabel',[]);
                            set(ax14,'xticklabel',[]);
                            eval(['ylabel(ax' num2str(12+n+(j==2)*1) ','' Velocity [km/s]'')']);
                            xlabel(ax13,'Time [ms]');
                            mBox2=uicontrol(h11,'Style','Text');% Make the title
                            set(mBox2,'String',['Raw Data and Reconstruction, Impact: ' num2str(dat(1).impacts(plotChan)) '[cm]']); 
                            set(mBox2,'fontweight','bold'); set(mBox2,'fontsize',13);
                            set(mBox2,'Position', [ 80,960,450,20]);set(mBox2,'BackgroundColor',[1 1 1]);
                            % Explainatory plot stuff
                            if j==2
                                plot(ax16,[reconTime(1)-.025,reconTime(end)+.025],ones(1,2).*param(plotChan,2+5*(i-1),n),'k--')
                                plot(ax16,ones(1,2)*(28/14500 +.25/14500 -dataPhase(plotChan,i)/(2*pi*14500))*1e3,[1,13],'k--')
                                title(ax16,['Raw Data and Fit, Shot: ' num2str(in(n).shot) ', Impact: ' num2str(dat(1).impacts(plotChan)) ' [cm]']);
                                %legend(ax16,{['Upper Fiber, RMS Error: ' num2str(RMS(1,plotChan,3,10)) ' [km/s]'];['Lower Fiber, RMS Error: ' num2str(RMS(2,plotChan,3)) ' [km/s]']},'Location','SouthWest','fontsize',10)
                                xlabel(ax16,'Time [ms]');
                            end
                                
                        
                        end
                        
                        % Explainatory reconsturction
                    else
                        plotChan=14;
                        eval(['plot(ax' num2str(12+n) ',dat(1).time.*in(n).timeScale,dat(in(n).line).vel(:,plotChan),''color'',partialColor{n,1},''linewidth'',lnwdth)']);
                        eval(['plot(ax' num2str(12+n) ',dat(1).time.*in(n).timeScale,param(plotChan,2,n)+param(plotChan,3,n)*sin(2*pi*dat(1).time.*1e-3.*in(n).timeScale.*14500 + dataPhase(plotChan)),''color'',partialColor{n,1},''linewidth'',lnwdth,''linestyle'',''--'')']);
                        set(ax15,'xticklabel',[]);
                        set(ax14,'xticklabel',[]);
                        eval(['title(ax' num2str(12+n) ', ''Shot: ' num2str(in(n).shot) ''')']);
                        eval(['ylabel(ax' num2str(12+n) ','' Velocity [km/s]'')']);
                        xlabel(ax13,'Time [ms]');
                        mBox2=uicontrol(h11,'Style','Text');% Make the title
                        set(mBox2,'String',['Raw Data and Reconstruction, Impact: ' num2str(dat(1).impacts(plotChan)) '[cm]']); 
                        set(mBox2,'fontweight','bold'); set(mBox2,'fontsize',13);
                        set(mBox2,'Position', [ 80,960,450,20]);set(mBox2,'BackgroundColor',[1 1 1]);
                    end
                    linkaxes([ax10,ax11,ax12],'x');
                    linkaxes([ax13,ax14,ax15],'x');
                    end
                    
                    %% HARDCODED DATA
                    if in(n).line==1 && strcmp(dat(1).title,'Shot 151217026'); % FOR 151217026 LINE 1
                       % dataPhase(11-1,1) = 6.14; % real weak signal, this is necessary because reasons
                        %dataPhase(21,1)= 3.5; % When signal is this weak, make it match the other fiber
                    elseif in(n).line==3 && strcmp(dat(1).title,'Shot 151217026')
                        if i==1
                            %dataPhase(:,1) = dataPhase(:,1)+2*pi;
                          %  dataPhase(2-1,1) = 2*pi;
                          %  dataPhase(4-1,1) = 5.8;
                        end
                    elseif in(n).line==2 && strcmp(dat(1).title,'Shot 151217026')
                        if i==2
                         %   dataPhase(7-1,2)=5;
                        end
                    elseif in(n).line==1 && strcmp(dat(1).title,'Shot 151217021')
                       % dataPhase(2,1)=0;
                        dataPhase(1,1)=3.2;
                    elseif (n==2 || n==3) && strcmp(dat(1).title,'Shot 160728011')
                        dataPhase(:,1) = dataPhase(:,1)+2*pi;
                    end
                    %dataPhase(:,i) = mod(dataPhase(:,i),5.43); % Minimum resolution, 151217026
                    
                  
                    
                    % half of integral of one half period is radius of motion
                    dataDispl(:,i) = param(:,3+5*(i-1),n).*(1./(param(:,5+5*(i-1),n)*2*pi)) .*1e5;

                else % if no FFT
                    i
                     % find index corresponding to time bounds
                    nTimeLim(1) = find(dat(1).time.*in(n).timeScale >= timebound(1), 1);
                    nTimeLim(2) = find(dat(1).time.*in(n).timeScale <= timebound(end), 1, 'last');
            %         nTimeLim(2) = nTimeLim(2) - 1; % the above command goes one too far

                    % Calculate all data
                    for m = 1:size(data, 2);
                        assignin('base','data',data);
                        selection = data((nTimeLim(1):nTimeLim(2))+(i-1)*size(data,1)/2, m);
                        dataAvg(m,i) = mean(selection(~isnan(selection)));
                        dataStd(m,i) = std(selection(~isnan(selection)));
                        dataDispl(m,i) = dataStd(m,i).*(1/14500)./(2*pi) .*1e5;

                    end
                end

                 % Plot Data
                 hold on;
                if Analysis==2 % plot the individual velocities
                    %for j=1:length(dataAvg(:,i));try dataAvg(j.*(pRel(j,i)<.5),i)=NaN;end;end
                    %t2(n) = errorbar(ax7,dat(1).impacts(1:size(data,2)), dataAvg(:,i), dataStd(:,i), 'color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
                elseif Analysis==1
                    t2(n) = errorbar(dat(1).impacts(1:size(data,2)), dataAvg(:,i), dataStd(:,i), 'color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
                end
                % Plot Displacement
                if in(n).fftPlot
                    
                    if Analysis==1
                        figure(h4);
                        hold on;
                        plot(dat(1).impacts(1:size(data,2)),dataDispl(:,i),'color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
                        ylabel('Average Radius of Displacement [cm]');set(gca,'ylim',[0,8]);
                        figure(h2);
                    elseif Analysis==2
                        for j=1:length(dataDispl(:,i));try dataDispl(j.*(pRel(j,i)<CutPow),i)=NaN;end;end
                        if ~plotError
                            saveDat(n).Displacement(:,i) = dataDispl(:,i);
                            plot(ax6,dat(1).impacts(1:size(data,2)),dataDispl(:,i),'-*','color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
                        else
                            error = sqrt( (nanmean(dat(in(n).line).velU(:,doubleplot(i,:)))./(2*pi*14500) .*1e5).^2 + squeeze(RMS(i,:,n,100)).^2) ; % convert error to centemeters displacement
                            if any(dispSupress(n,:,i))
                                dataDispl(dispSupress(n,:,i),i)=NaN;
                            end
                            saveDat(n).Displacement(:,i) = dataDispl(:,i);
                            saveDat(n).DisplacementError(:,i) = error;
                            errorbar(ax6,dat(1).impacts(1:size(data,2)),dataDispl(:,i),error,'-*','color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
                        end
                         ylabel(ax6,'[cm]');set(ax6,'ylim',[0,8]);
                        set(ax6,'xlim',xlim);
                        if ~includeTemp
                            xlabel(ax6,'Impacts [cm]');
                        else
                            xlabel(ax17,'Impacts [cm]');
                        end
                            
                    end
                end
            end
            
            %%%%%%%%%%%%%% Plotting Analysis %%%%%%%%%%%%%%%%%%%%%%%%
            if in(n).fftPlot & in(n).doubleplot==1 % This one is usually executed
                % plot Flow Profile
                for i = 1:size(data, 2)
                    cycle1 = (param(i,3,n)*sin(param(i,4,n)+(2*pi).*(0:(1/100):1))+param(i,2,n));
                    cycle2 = (param(i,8,n)*sin(param(i,9,n)+(2*pi).*(0:(1/100):1))+param(i,7,n));
                    maxDispl(i) = max(cycle1-cycle2);
                    minDispl(i) = min(cycle1-cycle2);
                end
                
                L=-(dataAvg(:,1)-dataAvg(:,2))+minDispl'; % This may be wrong. We might just want max/min Displ
                U=-(dataAvg(:,1)-dataAvg(:,2))-maxDispl';
                if Analysis==1
                    figure(h5)
                    errorbar(dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2))./2,L,U,'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style);
                    set(gca,'ylim',[-2,20]);
                    
                    figure(h3);
                elseif Analysis==2 
                    %for k=1:2;for j=1:length(dataAvg(:,k));try dataAvg(j.*(pRel(j,k)<.5),k)=NaN;end;end;end
%                      t3(n)=errorbar(ax7,dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2)),L,U,'color',['k'],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1},...
%                           'MarkerEdgeColor',[in(n).color{1}]);
                         if plotType==1
                             if any(flowSupress(n,:))
                                 dataAvg(flowSupress(n,:),1)=NaN;
                             end
                                 
                             if ~plotError
                                 saveDat(n).Flow = -(dataAvg(:,1)-dataAvg(:,2))./2;
                                 t3(n)=plot(ax7,dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2))./2,'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1},...
                             'MarkerEdgeColor',[in(n).color{1}]);
                             else
                                 error = sqrt( mean(dat(in(n).line).velU(:,doubleplot(1,:))).^2 + mean(dat(in(n).line).velU(:,doubleplot(2,:))).^2);
                                 saveDat(n).Flow = -(dataAvg(:,1)-dataAvg(:,2))./2;
                                 saveDat(n).FlowError = error;
                                 t3(n)=errorbar(ax7,dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2))./2,error,'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1},...
                        	'MarkerEdgeColor',[in(n).color{1}]);
                             end
                         elseif plotType==2
                             t3(n)=plot(ax7,dat(1).impacts(1:size(data,2)),dataAvg(:,1),'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1},...
                         'MarkerEdgeColor',[in(n).color{1}]);
                             t3(n)=plot(ax7,dat(1).impacts(1:size(data,2)),dataAvg(:,2),'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{2},...
                         'MarkerEdgeColor',[in(n).color{1}]);
                         end
                         
                       
                    set(ax7,'ylim',[-2,20]);
                end
                % Plot Phases
                hold on;
                if max(max(dataPhase))>2*pi
                    disp(['dataPhase > 2Pi, shifting by: ' num2str(-(max(max(dataPhase))-2*pi))]);
                    dataPhase=dataPhase-(max(max(dataPhase))-2*pi);
                end
                % can add 2Pi and maintain relationship with A injector
                if mean(mean(dataPhase)) <=-pi
                    dataPhase= dataPhase+2*pi;
                    disp('Shifting up by 2pi');
                end
                if mean(mean(dataPhase)) >=pi
                    dataPhase= dataPhase-2*pi;
                    disp('Shifting down by 2pi');
                end
                
                %% HARDCODING PHASE CHANGES, 2PI LEGAL
                if (n==3) && in(1).shot==160728011 && plotType ==1
                        dataPhase(:,1) = dataPhase(:,1)+2*pi;
                end
                if in(n).shot == 160728011 && in(n).line==2&& plotType ==1
                    dataPhase(11:17,1) = dataPhase(11:17,1)+2*pi;
                   dataPhase(3,1) = dataPhase(3,1)+2*pi;
                end
                
                if in(n).shot == 160728011 && in(n).line==1&& plotType ==1
                    dataPhase(12:16,1) = dataPhase(12:16,1)+2*pi;
                end
                if in(n).shot == 160728013 && in(n).line==2 && plotType ==1
                    dataPhase(1:10,1) = dataPhase(1:10,1)-2*pi;
                    dataPhase(:,1) = dataPhase(:,1)+2*pi;
                end
                if in(n).shot == 160728013 && in(n).line==1 && plotType ==1
                    dataPhase(11,1) = dataPhase(11,1)-2*pi;
                    dataPhase(12,1) = dataPhase(12,1)+2*pi;
                    dataPhase(:,1) = dataPhase(:,1)+2*pi;
                    %dataPhase(6,2) = dataPhase(6,2)+2*
                end
                if in(n).shot == 160728013 && in(n).line==2 && plotType ==2
                    dataPhase=dataPhase;
                end
                if in(n).shot == 160728011 && in(n).line==2 && plotType ==2
                    dataPhase=dataPhase;
                end
                if in(n).shot == 160728012 && in(n).line==2 && plotType ==2
                    dataPhase=dataPhase+2*pi;
                    dataPhase([3:8,15],2) = dataPhase([3:8,15],2)+2*pi;
                end
                 if in(n).shot == 160728013 && in(n).line==1&& plotType ==1
                    dataPhase(1:8,1) = dataPhase(1:8,1)-2*pi;
                    dataPhase(12,1) = dataPhase(12,1)-2*pi;
                end
                if in(n).shot == 160728021 && in(n).line==1
                    dataPhase = dataPhase+pi;
                    dataPhase(:,1) = dataPhase(:,1)-2*pi;
                end
                if in(n).shot == 160728024 && in(n).line==1
                    dataPhase = dataPhase-pi;

                end
                if in(n).shot == 160728024 && in(n).line==2
                    dataPhase(:,2) = dataPhase(:,2)-2*pi;

                end
                if in(n).shot == 160728024 && in(n).line==2
                    dataPhase = dataPhase+2*pi;
                end
                if in(n).shot == 160728012 && in(n).line==2
                    %dataPhase(:,2) = dataPhase(:,2)+2*pi;
                end
                if in(n).shot == 8160609010
                    dataPhase(1:6,2) = dataPhase(1:6,2)-2*pi;
                    dataPhase(24:33,2) = dataPhase(24:33,2)-2*pi;
                    dataPhase(24:36,1) = dataPhase(24:36,1)-2*pi;
                end
                if in(n).shot == 160525016
                    dataPhase(6:7,1) = dataPhase(6:7,1) -2*pi;
                end
                if in(n).shot == 160525017
                    dataPhase(1:6,1) = dataPhase(1:6,1) -2*pi;
                    dataPhase(1:2,2) = dataPhase(1:2,2) -2*pi;
                end
                if in(n).shot == 160525018
                   dataPhase([5,6],2) = dataPhase([5,6],2) -2*pi;
                end
                if in(n).shot == 160525019
                   dataPhase(5:6,2) = dataPhase(5:6,2) -2*pi;
                end
                if any([in.shot] == 160728021)
                    injParam(3)=injParam(3)-2*pi;
                end
                
                % Nimrod HITSI
                if in(n).shot == 8129499
                   dataPhase(1:19,:) = dataPhase(1:19,:) +2*pi;
                   %dataPhase(:,2) = dataPhase(:,2) -2*pi;
                   dataPhase(16:19,2) = dataPhase(16:19,2) -2*pi;
                end
                if in(n).shot==151217024
                    dataPhase = dataPhase-2*pi;
                    dataPhase(15,2)=dataPhase(15,2)-2*pi;
                end
                if in(n).shot==151217025
                    dataPhase(15:17,2) = dataPhase(15:17,2)-2*pi;
                end
                % Low Perf HitSI3
                if in(n).shot == 151217026
                    dataPhase(4:12,1)=dataPhase(4:12,1)-2*pi;
                    dataPhase(7:17,2)=dataPhase(7:17,2)-2*pi;
                end
                
                %% Flip upper array, look at when toroidal drive happens
                if  driveDirection  == 1
                        dataPhase(:,1) = dataPhase(:,1)-pi;
                        display('WARNING: POSITIVE FIBER DROPPED -PI, DRIVE DIRECTION PHASE');
                end
                dataPhase = dataPhase +in(n).phaseShift;
                
                for i=1:size(doubleplot,1)
                    if Analysis==1
                        plot(dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,'-*','color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
                        ylabel('Phase [deg]');
                    elseif Analysis==2 && in(1).shot > 8129499
                        for j=1:length(dataPhase(:,i));try dataPhase(j.*(pRel(j,i)<CutPow),i)=NaN;end;end
%                         ax8Data = dataPhase(:,i);
%                         if plotType==1
%                             xindex=dat(1).impacts(1:size(data,2));
%                         elseif plotType==2
%                             xindex=dat(1).impacts(1*(i==1) +(size(data,2)+1)*(i==2):size(data,2)*(i==1)+(end)*(i==2));
%                         end
                        if any(phaseSupress(n,:,i))
                            dataPhase(phaseSupress(n,:,i),i)=NaN;
                        end
                        if ~plotError
                            saveDat(n).Phase(:,i) = dataPhase(:,i).*180./pi;
                            phaseH(:,i,n)=plot(ax8,dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,'-*','color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
                        else
                            saveDat(n).Phase(:,i) = dataPhase(:,i).*180./pi;
                            saveDat(n).PhaseError(:,i) = SigDev(n,:,i).*180./pi;
                            phaseH(:,i,n)=errorbar(ax8,dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,SigDev(n,:,i).*180./pi,'-*','color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
                        end                            
                        ylabel(ax8,'[deg]');set(ax8,'ylim',[-400,400]);set(ax8,'xticklabel',[]);
                        set(ax8,'xlim',xlim);
                        if n==1
                            saveDat(n).injPhase = injParam(:,3).*180./pi;
                            saveDat(n).injPhase(1) = saveDat(n).injPhase(1)-20;
                        plot(ax8,[-60 60]',[injParam(1,3) injParam(1,3)]'.*180./pi-20,'--','LineWidth',lnwdth,'color',[0    0.4470    0.7410]);
                        plot(ax8,[-60 60]',[injParam(2,3) injParam(2,3)]'.*180./pi,'--','LineWidth',lnwdth,'color',[0.8500    0.3250    0.0980]);
                        plot(ax8,[-60 60]',[injParam(3,3) injParam(3,3)]'.*180./pi,'--','LineWidth',lnwdth,'color',[0.9290    0.6940    0.1250]);
                        if in(1).shot >= 160728011 || in(1).shot == 151217024
                            plot(ax8,[-60 60]',[injParam(3,3) injParam(3,3)]'.*180./pi + 360,'--','LineWidth',lnwdth,'color',[0.9290    0.6940    0.1250]);
                        end
                        end
                        
                    elseif Analysis==2 && in(1).shot <= 8129499
                        for j=1:length(dataPhase(:,1));try dataPhase(j.*(pRel(j,i)<CutPow),1)=NaN;end;end
                        phaseH(1,1,n)=plot(ax8,dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
                        ylabel(ax8,'[deg]');set(ax8,'ylim',[0,400]);set(ax8,'xticklabel',[]);
                        set(ax8,'xlim',xlim);
                        if n==1
                        plot(ax8,[0 60]',[injParam(1,3) injParam(1,3)]'.*180./pi,'--','LineWidth',lnwdth,'color',[    0    0.4470    0.7410]);
                        plot(ax8,[0 60]',[injParam(2,3) injParam(2,3)]'.*180./pi,'--','LineWidth',lnwdth,'color',[    0.8500    0.3250    0.0980]);
                        end

                    end
                    
                    % Plot temperature 
                    if includeTemp
                        if any(tempSupress(n,:,i))
                            dat(in(n).line).temp(:,doubleplot(i,tempSupress(n,:,i)))=NaN;
                        end
                        if ~plotError
                            saveDat(n).Temp(:,i) = mean(dat(in(n).line).temp(:,doubleplot(i,:)));
                            plot(ax17,dat(1).impacts(1:size(data,2)),mean(dat(in(n).line).temp(:,doubleplot(i,:))),'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
                        else
                            saveDat(n).Temp(:,i) = mean(dat(in(n).line).temp(:,doubleplot(i,:)));
                            saveDat(n).TempError(:,i) = mean(dat(in(n).line).tempU(:,doubleplot(i,:)));
                            errorbar(ax17,dat(1).impacts(1:size(data,2)),mean(dat(in(n).line).temp(:,doubleplot(i,:))),mean(dat(in(n).line).tempU(:,doubleplot(i,:))),'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
                        end
                        ylabel('[eV]');
                        set(ax6,'xticklabel',[]);
                    end
                end
                
                
                %  Plot Fft Spectrum
                plot(ax9,dat(1).impacts(1:size(data,2)),100*pRel(:,1),'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
                plot(ax9,dat(1).impacts(1:size(data,2)),100*pRel(:,2),'-*','color', in(n).color{2}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{2});
                 ylabel(ax9,'[%]');set(ax9,'ylim',[0,100]);
                 set(ax9,'xlim',xlim); title(ax9,'Inj. Mode % of Reconstruction');
                plot(ax9,xlim,[CutPow,CutPow].*100,'--k');
                if n==1;figure(h2);end
                
            %%%%%%%%%%%%%%%HITSI%%%%%%%%%%%%%%%%
            
            elseif ~isempty(in(n).fftPlot) && Analysis ==2 
                display('EXECUTING HITSI DATA')
                % Plot Displacement Profile
                %if in(n).shot > 130000
%                 t3(n)=errorbar(ax7,dat(1).impacts(1:size(data,2)),param(:,2,n),param(:,3,n),param(:,3,n),'color',['k'],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1},...
%                     'MarkerFaceColor',[in(n).color{1}]);
                % Check data supression
               if any(flowSupress(n,:))
                    param(flowSupress(n,:),3,n)=NaN;
                end
                if ~plotError
                     saveDat(n).Flow = param(:,2,n);
                     t3(n)=plot(ax7,dat(1).impacts(1:size(data,2)),param(:,2,n),'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1});%,...
                else
                    saveDat(n).Flow = param(:,2,n);
                     saveDat(n).FlowError = nanmean(dat(in(n).line).velU);
                     t3(n)=errorbar(ax7,dat(1).impacts(1:size(data,2)),param(:,2,n),nanmean(dat(in(n).line).velU),'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1});%,...
                end
                    %'MarkerEdgeColor',[in(n).color{1}]);
                %end
                set(ax7,'ylim',[-2,15]);
               
                % Plot Phases
                hold on;
                if max(max(dataPhase))>2*pi
                    dataPhase=dataPhase-(max(max(dataPhase))-2*pi);
                end
                for j=1:length(dataPhase(:,1));try dataPhase(j.*(pRel(j,i)<CutPow),1)=NaN;end;end
                %% HARDCODING
                if in(1).shot==129499 && n==2 
                    %dataPhase(8:end)=dataPhase(8:end)+2*pi;
                end
                if in(1).shot==129499 && n==1
                    %dataPhase(9:(end-3))=dataPhase(9:(end-3))+2*pi;
                end
                dataPhase = dataPhase + in(n).phaseShift;
               if any(phaseSupress(n,:))
                    dataPhase(phaseSupress(n,:),1)=NaN;
                end
                if ~plotError
                    saveDat(n).Phase(:,1) = dataPhase(:,1).*180./pi;
                    phaseH(1,1,n)=plot(ax8,dat(1).impacts(1:size(data,2)),dataPhase(:,1).*180./pi,'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
                else
                    saveDat(n).Phase(:,1) = dataPhase(:,1).*180./pi;
                    saveDat(n).PhaseError(:,1) = SigDev(n,:,1).*180./pi;
                    phaseH(1,1,n)=errorbar(ax8,dat(1).impacts(1:size(data,2)),dataPhase(:,1).*180./pi,squeeze(SigDev(n,:)).*180./pi,'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
                end
                ylabel(ax8,'[deg]');set(ax8,'ylim',[0,400]);set(ax8,'xticklabel',[]);
                set(ax8,'xlim',xlim);
                % Plot injector Phase
                if n == 1
                saveDat(n).injPhase = injParam(:,3).*180./pi;
                plot(ax8,[0 60]',[injParam(1,3) injParam(1,3)]'.*180./pi,'--','LineWidth',lnwdth,'color',[    0    0.4470    0.7410]);
                plot(ax8,[0 60]',[injParam(2,3) injParam(2,3)]'.*180./pi,'--','LineWidth',lnwdth,'color',[    0.8500    0.3250    0.0980]);
                end
                % Plot Temperature
                if any(tempSupress(n,:))
                    dat(in(n).line).temp(:,tempSupress(n,:))=NaN;
                end
                if includeTemp
                    if ~plotError
                        saveDat(n).Temp = nanmean(dat(in(n).line).temp);
                        plot(ax17,dat(1).impacts(1:size(data,2)),nanmean(dat(in(n).line).temp),'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
                    else
                        saveDat(n).Temp = nanmean(dat(in(n).line).temp);
                        saveDat(n).TempError = nanmean(dat(in(n).line).tempU);
                        errorbar(ax17,dat(1).impacts(1:size(data,2)),mean(dat(in(n).line).temp),nanmean(dat(in(n).line).tempU),'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
                    end
                    ylabel(ax17,'[eV]');
                    set(ax6,'xticklabel',[]);
                end
                 %  Plot Fft Spectrum
                plot(ax9,dat(1).impacts(1:size(data,2)),100*pRel(:,1),'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
                plot(ax9,dat(1).impacts(1:size(data,2)),100*pRel(:,2),'-*','color', in(n).color{2}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{2});
                 ylabel(ax9,'[%]');set(ax9,'ylim',[0,100]);
                 set(ax9,'xlim',xlim); title(ax9,'Inj. Mode % of Reconstruction');
                plot(ax9,xlim,[CutPow,CutPow].*100,'--k');
                if n == 1;figure(h2);end
            end
            % Plot Flow Profile
            %plot(dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2)),'k','LineWidth', lnwdth, 'LineStyle', in(n).style);
            if Analysis==1 & in(n).doubleplot==1
                
                plot(dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2))./2,'color' ,'k','marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1});
                plot(xlim,[0,0],'--k')
                ylabel('Toroidal Flow [km/s]'); set(gca,'ylim',[-10,10]);
                xlabel('Impacts [cm]');
            elseif Analysis==2 %& in(n).doubleplot==1
               % plot(ax7,dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2)),'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style);
                plot(ax7,xlim,[0,0],'--k'); set(ax7,'xlim',xlim);
                ylabel(ax7,'[km/s]'); set(ax7,'ylim',[-6,6]);set(ax7,'xticklabel',[]);
                h2.delete;

                linkaxes([ax6 ax7 ax8 ],'x');
                if includeTemp
                    linkaxes([ax6 ax7 ax8 ax17 ],'x');
                end
            end

            figure(h) % make other current
        elseif plotType ==2
            % find cycle averaged temperature and fluxuation amplitude
            
        
        
        end
    end
    
    %% offset each line for plot 1
    for j = 1:size(data, 2)
       %PhaseVelocity.Velocity(j,:) = PhaseVelocity.Velocity(j,:)+ (j-1)*50;   
        data(:, j) = data(:, j) + (j-1) * offset;
        zeroline(:,j) = zeros(ceil(size(data,1)/2),1)+(j-1) * offset;
    end

    %% Plot Data
%     if in(n).error
%         time = ndgrid(dat(1).time, 1:size(data, 2));
%         t(n, :) = errorbar(ax, time, data, errorL, errorU, 'color', in(n).color, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
% %         for m = 1:size(t, 2)
% %             errorbar_tick(t(n, m), errWdth); % adjust errorbar width
% %         end
%     else
        if ~isempty(in(n).doubleplot) && ~plotTor % plot both fibers
            saveDat(n).UpperFiberData = data(1:length(dat(1).time),:);
            saveDat(n).LowerFiberData = data(length(dat(1).time)+1:end,:);
            saveDat(n).Impacts = dat(1).impacts(1:size(data,2));
            saveDat(n).Time =time(1:length(dat(1).time));
                if ~in(n).error
                t(n, :) = plot(ax, time(1:length(dat(1).time)), data(1:length(dat(1).time),:), ...
                    'color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
                t(n, :) = plot(ax, time(length(dat(1).time)+1:end)', data(length(dat(1).time)+1:end,:), ...
                    'color', in(n).color{2}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{2});
            
            else
                 t(n, :) = errorbar(ax, time(1:length(dat(1).time)), data(1:length(dat(1).time),:),...
                     dat(in(n).line).velU(:,doubleplot(1,:)),dat(in(n).line).velL(:,doubleplot(1,:)),...
                    'color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
                t(n, :) = errorbar(ax, time(length(dat(1).time)+1:end)', data(length(dat(1).time)+1:end,:), ...
                    dat(in(n).line).velU(:,doubleplot(1,:)),dat(in(n).line).velL(:,doubleplot(1,:)),...
                    'color', in(n).color{2}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{2});
            end 
            plot(ax, [time(1),time(length(dat(1).time))], zeroline([1,size(data,1)./2],:), ...
                '--k', 'LineWidth', .5, 'LineStyle', in(n).style{1});
        elseif  ~isempty(in(n).doubleplot) && plotTor % plot the toroidal flow
            t(n, :) = plot(ax, time(1:length(dat(1).time)), data(length(dat(1).time)+1:end,:)-data(1:length(dat(1).time),:) +zeroline, ...
                'color', 'k', 'LineWidth', lnwdth, 'LineStyle', in(n).style);
            plot(ax, [time(1),time(length(dat(1).time))], zeroline([1,size(data,1)./2],:), ...
                '--k', 'LineWidth', .5, 'LineStyle', in(n).style);
        else
             saveDat(n).UpperFiberData = data(1:length(dat(1).time),:);
             saveDat(n).Impacts = dat(1).impacts(1:size(data,2));
            saveDat(n).Time =time(1:length(dat(1).time));
            if ~in(n).error
                t(n, :) = plot(ax, time, data, 'color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
            else
                size(data)
                size(dat(in(n).line).velU)
                size(dat(in(n).line).velL)
                t(n, :) = errorbar(ax, time, data,dat(in(n).line).velU,dat(in(n).line).velL,...
                'color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
            end
        end
    end
    
% end
if ~isempty(in(1).fftPlot)
    % Make the Legend
    [leg,icon] = legend(ax8,phaseH(1,1,:),{in(1:1:end).legend});
    for i = 0: length(in)-1
        try %(~plotError || in(1).shot < 829499)
            set(icon(length(in)+1+2*i),'Color',in(i+1).color{1});
            set(icon(length(in)+2+2*i),'Color',in(i+1).color{1});
            set(icon,'LineStyle','-');
        catch
            set(icon(length(in)+1+i).Children.Children,'Color',in(i+1).color{1});
            set(icon(i+(floor(length(icon))/2+1)).Children.Children,'LineStyle','-');
        end
    end
    


%     % Annotation
%     figure(h6)
%     a=annotation('textbox',[.2,.6,.3,.3],'String',{'Upper Array: Solid','Lower Array: Dashed'},'FitBoxToText','on');
%     %a.Position=[.175,.84,.175,.075];
%     a.Position=[.175,.81,.175,.11];
%     
%     a1=annotation('textbox',[.2,.6,.3,.3],'String',{'1/2 Difference Between Upper and Lower Arrays'},'FitBoxToText','on');
%     %a1.Position=[.175,.54,.2,.075];
%     a1.Position=[.175,.50,.2,.11];
%     
%     a2=annotation('textbox',[.2,.6,.3,.3],'String',{'Amplitude of Sine Fit'},'FitBoxToText','on');
%     %a2.Position=[.175,.27,.175,.045];
%     a2.Position=[.175,.25,.175,.06];
%     
%     a4=annotation('textbox',[.2,.6,.3,.3],'String',{'Temporal Phase Offset for Sine Fit'},'FitBoxToText','on');
%     %a4.Position=[.175,.75,.175,.06];
%     a4.Position=[.175,.72,.175,.08];

end

figure(h); % make first figure current

%% Impact Parameter Labels on Right
for n = 1:size(data, 2)
    y = offset * (n-1) + 0.1 * offset;
    text(timebound(end) + 0.02 * (timebound(end) - timebound(1)), y, num2str(dat(1).impacts(n), 2), 'fontsize', fntsz);
    plot(ax, [dat(1).time(1), dat(1).time(end)], zeros(2) + (n-1) * offset, '-', 'color', 'k');
end

% %% Legend
% legendText = cell(1, length(in)); % initialize
% legendHands = zeros(1, length(in)); % initialize
% for n = 1:length(in)
%     legendText{n} = in(n).legend;
%     legendHands(n) = t(n, 1);
% end
% legend(ax, legendHands, 'Location', 'NorthEastOutside',...
%      legendText, 'fontsize', fntsz);

 %% Misc. Figure Properties
if ~compactCurrents
    xlabel('Time [ms]','fontsize', fntsz);
else
    set(gca,'xticklabel',[]);
end
ylabel(sidebar, 'fontsize', fntsz);

title([in(1).AnalysisTitle  ':' titles], 'fontsize', fntsz);
%set(ax, 'XLim', timebound);
set(gca, 'LineWidth', lnwdth);
set(gca, 'fontsize', fntsz);
box on;
grid on;
if or(plotType == 2, plotType == 3) % temperature or Intensity
    yLowerLim = 0;
    yUpperLim = offset * size(data, 2);
else
    yLowerLim = -0.5 * offset;
    yUpperLim = offset * size(data, 2) - 0.5 * offset;
end
set(gca, 'YLim', [yLowerLim, yUpperLim]);
set(gca, 'YTick', []);

nt = text(timebound(end) + 0.125 * (timebound(end) - timebound(1)), (offset * size(data, 2))/2, 'R  [cm]', 'fontsize', fntsz);
set(nt, 'rotation', -90)

%% Figure Properties for Averages
if plotAverages
    figure(h2) % make current
    
    title(['Mean ' titles ' and Fluctuations'], 'fontsize', fntsz);
    set(ax2, 'XLim', [dat(1).impacts(1) - 1, dat(1).impacts(end) + 1]);
    set(gca, 'LineWidth', lnwdth);
    set(gca, 'fontsize', fntsz);
    box on;
    grid on;
    xlabel('Major Radius [cm]');
    ylabel(units);
    
    % Legend
    legendText = cell(1, length(in)); % initialize
    legendHands = zeros(1, length(in)); % initialize
    for n = 1:length(in)
        legendText{n} = in(n).legend;
        legendHands(n) = t2(n);
    end
    legend(ax2, legendHands, 'Location', 'NorthEastOutside',...
        legendText, 'fontsize', fntsz);
end

%% Plot Currents
if plotCurrents
    % plot currents
    if ~compactCurrents
        j = figure('Visible', 'on', 'Name', 'MULTIPLOT-Currents', 'Position',...
            [5, 35, figureWidth, analysisHeight], 'Color', [1 1 1]);
        ax2 = axes('Parent', j, 'Position', [0.15, 0.08, 0.8, 0.15], 'FontSize', fntsz+1); % currents
    else
        ax2 = axes('Parent', h, 'Position', [.075, .08, .8, .15], 'FontSize', fntsz+1);
    end
    if length(in) > 1
        load(['dat' num2str(in(1).shot) '10.mat']);
    end
    if dat(1).shotRef >999999
        plot(ax2, dat(1).iinjaTime.*in(1).injTimeScale, ...
            dat(1).iinja.*in(1).injScale,'LineWidth', lnwdth','color',[0    0.4470    0.7410]); 
        hold on;
        plot(ax2, dat(1).iinjbTime.*in(1).injTimeScale, ...
            dat(1).iinjb.*in(1).injScale, 'color',[0.8500    0.3250    0.0980], 'LineWidth', lnwdth);
        hold on;
        plot(ax2, dat(1).iinjcTime.*in(1).injTimeScale, ...
            dat(1).iinjc.*in(1).injScale, 'color',[0.9290    0.6940    0.1250], 'LineWidth', lnwdth);
        hold on;
    else
         plot(ax2, dat(1).iinjxTime.*in(1).injTimeScale, ...
            dat(1).iinjx.*in(1).injScale,'LineWidth', lnwdth,'color',[0    0.4470    0.7410]); 
        hold on;
        plot(ax2, dat(1).iinjyTime.*in(1).injTimeScale, ...
            dat(1).iinjy.*in(1).injScale, 'color',[0.8500    0.3250    0.0980], 'LineWidth', lnwdth);
        hold on
    end
    for i = 1:length(in)
        plot(ax2, dat(1).ItorTime.*in(1).injTimeScale, ...
            Itor(:,i).*in(1).injScale, 'color', in(i).color{1}, 'LineWidth', lnwdth);
    end
    set(gca, 'LineWidth', lnwdth);
    xlabel('Time [ms]');
    ylabel('I_{INJ} [kA]');
    linkaxes([ax,ax2],'x');
    set(ax2, 'XLim', timebound);
    grid on;
    
end
cd(['T:\IDS\Analysis Repository\']);
% if in(1).shot == 151217023 || in(1).shot == 151217024 || in(1).shot == 151217025 || in(1).shot == 151217026
%     cd(['T:\IDS\Analysis Repository\' num2str(in(1).shot) ]);
% else
%     cd(['T:\IDS\Analysis Repository\Alternate Phasings\60D2\']);
% end
%saving
if saving
    pause(1.5);% necessary to get figure size correct
    saveas(h, ['Lines L' num2str(in(1).line) titles num2str(in(1).shot)], 'fig');
    saveas(h6, ['Analysis L' num2str(in(1).line) num2str(in(1).shot)], 'fig');
    saveas(h, ['Lines L' num2str(in(1).line) titles num2str(in(1).shot)], 'bmp');
    saveas(h6, ['Analysis L' num2str(in(1).line) num2str(in(1).shot)], 'bmp');
    if in(1).fftPlot
        saveas(h7, ['Spectrum L' num2str(in(1).line) num2str(in(1).shot)], 'fig');
        saveas(h7, ['Spectrum L' num2str(in(1).line) num2str(in(1).shot)], 'bmp');
    end
end

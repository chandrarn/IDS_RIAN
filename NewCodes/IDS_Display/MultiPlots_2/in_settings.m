%% Input Settings for Multiplots
% Returns the 'in' structure, which includes the shots, formatting stuff,
% and settings for plotting the FFT, etc
% Also returns the range of channels in each array to plot, and the time
% range to plot over

function [in,timebound,chan_range,xlim,supress] = in_settings(n)
lines = {'O II', 'C III', 'O II','C III'};

switch n
    case 1
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
        % 
        in(1).shot = 160728013;
        in(1).line=2;
        in(1).color = {'r';'r'};%[66, 188, 244]./255};
        in(1).style = {'-','--'};
        in(1).legend = [num2str(in(1).shot) ': +65.6kA'];
        in(1).phaseShift = -pi/2;
        in(1).error = 0; % 1 / 0 for errorbars
        in(1).velShift = -5; % SHIFT VELOCITY
        in(1).intScale = 2; % scale factor for intensity
        in(1).timeShift = 0; % ms, shift time base
        in(1).timeScale = 1e-3; % scale timebase to put into ms
        in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
        in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
        in(1).doubleplot = [1];[1:23; 24,26:47]; % plot coorespoinding impacts
        in(1).fftPlot = [1]; % FFT of signal, n frequencies
        in(1).AnalysisTitle=[ 'HIT-SI3: 0-120-240 Phasing, Gain 2.9, ' lines{in(1).line}];
        in(1).phaseShift = 2*pi -pi/2 ;% The +pi is to overlay with the positive shot
        in(2)=in(1);
        in(2).phaseShift = 2*pi -pi/2-pi;
        in(2).shot = 160728012;
        in(2).color = {'b';'b'};%[182, 244, 66]./255};
        in(2).legend = [num2str(in(2).shot) ': -52.3kA'];
        %Supression of Dissenting data:
        % For C III
        if in(1).line==2
        % Format: Shot - Index from Outboard - Array
        phaseSupress(1,:,2) = [4];
        phaseSupress(1,:,1)=1;
        phaseSupress(2,:,2)=0;
        flowSupress(1,:)=[1];
        flowSupress(2,:)=[1];
        dispSupress(1,:,2)=[1,1,1,1,1];
        dispSupress(2,:,2)=[1,1,1,1,4];
        dispSupress(2,:,1)=[1,1,1,1,1];
        tempSupress(1,:,1)=[2,2,2,2];
        %tempSupress(1,:,2)=[1:2,5];
        tempSupress(1,:,2)=[1,2,2,2];
        %tempSupress(2,:,2)=[1,3,5];
        tempSupress(2,:,2)=[1,1,1,1];
        tempSupress(2,:,1)=[1,1,19,1];
        else
        phaseSupress(1,:,1)=[18,19];
        phaseSupress(2,:,2)=[0,0];
        flowSupress(1,:)=[1,2];
        flowSupress(2,:)=[1,1];
        dispSupress(1,:,2)=[1,1,1,2,3];
        dispSupress(2,:,2)=[1,1,3,2,4];
        dispSupress(2,:,1)=[1,1,1,1,1];
        tempSupress(1,:,1)=[2,3,4,5];
        %tempSupress(1,:,2)=[1:2,5];
        tempSupress(1,:,2)=[1,2,2,3];
        %tempSupress(2,:,2)=[1,3,5];
        tempSupress(2,:,2)=[1,2,1,3];
        tempSupress(2,:,1)=[1,1,19,3];
        end

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


    case 2
        %% Low Performance
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

    case 3
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


    case 4
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

    case 5
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

    case 6
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


    case 7
        % % 129499 
        % % note: need to change .*1e-6 to .*1e-3 in sinefit
        % 
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
        in(1).fftPlot = [1]; % FFT of signal, n frequencies
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
        in(2).line=2;%1 is OII 2 is CIII
        in(2).shot = 129496;
        in(2).color = {'b';'b'};%[182, 244, 66]./255};
        in(2).legend = [num2str(in(2).shot) ': -76kA'];
        % Supression of Dissenting data:
        % for O II
        phaseSupress(1:2,:)=22;
        flowSupress(1:2,:)=[NaN,NaN;22,22];
        dispSupress(1:2,:)=22;
        tempSupress(1:2,:)=[21,22;21,22];
        % 
        % in(1).shot = 129499;%150625998;
        % in(1).line = 3; % line # NB: 1 is C III, 2 is O II, 3 is C III !
        % in(1).legend = [num2str(in(1).shot) ': +90kA'];
        % in(1).color = {'r';'r'};%[225,105,0]./255};
        % in(1).style = {'-','--'};
        % in(1).error = 0; % 1 / 0 for errorbars
        % in(1).velShift = 0; % SHIFT VELOCITY
        % in(1).intScale = 1; % scale factor for intensity
        % in(1).timeShift = 0; % ms, shift time base
        % in(1).timeScale = 1;1e-3; % scale timebase to put into ms
        % in(1).injTimeScale = 1;1e-3; % scale the injector time to ms
        % in(1).injScale = 1e0; 1e-3; % scale the inj current into kA
        % in(1).doubleplot = [];[1:23; 24,26:47]; % plot coorespoinding impacts
        % in(1).fftPlot = []; % FFT of signal, n frequencies
        % in(1).AnalysisTitle=['HIT-SI: 0-90 Phasing, ' lines{in(1).line+1}]; 
        % in(1).phaseShift =-pi/2;
        % in(1).error=0;
        % % 129450, 129451
        % % in(2) = in(1);
        % % in(2).shot=129450;
        % % in(2).line=1;
        % % in(2).legend = [num2str(in(2).shot) ': +79kA'];
        % % in(2).color = {'b';'b'};%[66, 188, 244]./255};
        % in(2)=in(1);
        % in(2).line=2
        % in(2).shot = 129496;
        % in(2).color = {'b';'b'};%[182, 244, 66]./255};
        % in(2).legend = [num2str(in(2).shot) ': -76kA'];
        % % Supression of Dissenting data:
        % phaseSupress(1:2,:)=22;
        % flowSupress(2,:)=22;
        % dispSupress(1:2,:)=22;
        % tempSupress(1:2,:)=[21,22;21,22];


    case 8
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

    case 9
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

        
    case 10
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

    case 11
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

end

supress.Phase = phaseSupress;
supress.Temp = tempSupress;
supress.Flow = flowSupress;
supress.Disp = dispSupress;



timebound = [1.1, 2.0]; % [ms]
% timebound = [1.2 2.3]; % [ms]

% Set these: ( add :2: if you only want every other line
% chan_range = 50:58;
% defaults
chan_ranget = [7:26];
chan_rangep = [43:62];
xlim = [0,50];


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

end
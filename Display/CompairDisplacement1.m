%% Displacement scatter plot, as a function of current amplification,
%  Compairing 14.5, 53.5, 63.5 kHz, with error bars
%  - Pull in data, average velocity, multiply by interval.
%  - Get std, get phase, get current amp. Plot.

clear all; %close all;
import MDSplus.*
HitTree = Connection('landau.hit');

% Frequency or Current Gain Plot?
Gain = 1;
% Pos/Neg torroidal current?
Pos = 1;

%Shots and parameters: These are hardcoded in based on the EPR poster
if ~Pos % negative torroidal current
    % 53&68 shots to use ( will find 14)
    shots(1).shots = [128585,128594,128592];% THE ORDER IS IMPORTANT:
    shots(2).shots = [129214,129213,129215];% IT MAKES PHASE-COLORS ~MATCH
    %shots(3).shots = [Velocity(:,2),Velocity(:,4),Velocity(:,10)]; % HARD CODED ORDER: MANUAL SELECTION
    %Load 14.5 velocity separated by phase:
    load('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data\Temp\Phase129817.mat');
    Velocity = PhaseVelocity.Velocity;
    Std = PhaseVelocity.Std;
    Std = [Std(:,2) Std(:,4) Std(:,10)];
    shots(3).shots = [Velocity(:,2) Velocity(:,4) Velocity(:,10)];
   
    %Time bounds ( Don't change 14.5 bounds: the averageing has already
    %occcured with the current time limits )
    TBase(1).time = [1.26, 1.95; 1.44, 1.99; 1.14, 1.53];
    TBase(2).time = [1.44, 1.84; 1.17, 1.77; 1.4,  1.83];
    TBase(3).time = [1.39, 1.73];

    % find the zero point, to counteract the offset:
    % WARNING: THIS MIGHT NOT BE SCIENCE
    % Calculated by averaging the most extreme datapoints.
    % For 53, this may be off by 25%
    Offset68 = 10.8050466248759; % average of 68(2) and 68(3)
    Offset53 = -3.27336872971015;% value for 53(3) ( This makes some sense, based on the relative phases and displacements )
    Offset14 = 1.97318465877948; % PRECISE MEAN DIFFICULT TO DETERMINE
        
else % Positive current
    % NOTE: no Pos 68 exists, using neg 68 
    % 53&68 shots to use ( will find 14)
    shots(1).shots = [128581,128593,128587,128586,128580];% THE ORDER IS IMPORTANT:
    shots(2).shots = [129214,129213,129215];% IT MAKES PHASE-COLORS ~MATCH
    %shots(3).shots = [Velocity(:,2),Velocity(:,4),Velocity(:,10)]; % HARD CODED ORDER: MANUAL SELECTION
    %Load 14.5 velocity separated by phase:
    load('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data\Temp\Phase129810.mat');
    Velocity = PhaseVelocity.Velocity;
    Std = PhaseVelocity.Std;
    Std = [Std(:,1) Std(:,5) Std(:,7) Std(:,8) Std(:,10)];
    shots(3).shots = [Velocity(:,1) Velocity(:,5) Velocity(:,7) Velocity(:,8) Velocity(:,10)];
   
    %Time bounds ( Don't change 14.5 bounds: the averageing has already
    %occcured with the current time limits )
    TBase(1).time = [1.19, 1.99; 1.39, 1.99; 1.19, 1.99; 1.5, 1.99; 1.45, 1.99];
    TBase(2).time = [1.44, 1.84; 1.17, 1.77; 1.4,  1.83];
    TBase(3).time = [1.47, 1.80];

    % find the zero point, to counteract the offset:
    % WARNING: THIS MIGHT NOT BE SCIENCE
    % Calculated by average shots ~ 180 out of phase:
    % assumption: displacement osculates around some zero point ( no net )
    Offset68 = 10.8050466248759; % average of 68(2) and 68(3)
    Offset53 = -4.59066437558201;% average of 53(2) and 53(3)
    Offset14 = -14.5127372746208;% Mean of dataset, std .55 ( they're all relatively close to the mean )
end

%Constants:
Offset = [Offset53 Offset68 Offset14];
%Frequencies
freq = [53.5, 68.5, 14.5];


torPlot = 0; % Poloidal plot
chan_ranget = [8:24]; % toroidal, mohawk port in midplane
% chan_ranget = [8:28]; % toroidal, mohawk port perp.
% chan_ranget = [8:27]; % toroidal, 71 degree port
% chan_ranget = [8:24]; % toroidal, axial port
% chan_ranget = 1:30; % NIMROD mohawk

% chan_rangep = [46:63]; % poloidal
% chan_rangep = [47:58]; % poloidal
chan_rangep = [50:58]; % poloidal, zoomed in on spheromak zone

if torPlot
    chan_range = chan_ranget;
else
    chan_range = chan_rangep;
end
    
% colors
ColorVector=hsv(10);

% markers
Marker = {'p';'s';'d'};

%TimeBase for CurrAmp
ts = -0.0001; % start time for new time base
%ts2 = ts - 0.0002; % start time for cropping data in TEST LATER
te = 0.005; % end time for new time base
%te2 = te + 0.0002; % end time for cropping data in
dtn = 0.4e-6; % dt for new time base (2.5 MHz)
tbase = ts: dtn: te; % new time base

%h1 = figure('Name','Displacement Comparison'); hold on;
figure;
ax = axes('FontSize',18,'LineWidth',3); hold on; grid on;

for i = 1:(2-Pos); % loop through frequencies. If positive, dont plot 68.5
    for j = 1: length(shots(i).shots)
        timeBounds = [TBase(i).time(j,1) TBase(i).time(j,2)]; % refresh time bounds
        % get the dat from the file
        load(['T:\IDS\Data Repository\TEMP\dat' num2str(shots(i).shots(j)) '10.mat']);
        
        %trim to the channel and time window we want
        while any(any(isnan(dat.vel)))
            dat = trimRandT(dat,chan_range,timeBounds);
            timeBounds(2)=timeBounds(2)-.01;
        end 
        clear t dt sig
        
        % Get Current Amplification
        HitTree.openTree('hitsi',shots(i).shots(j));% open shot
        % Get data
        [t.itorsm, dt.itorsm,sig.itorsm] = gen_data_in('sihi_smooth(\i_tor_spaavg)',HitTree);
        [t.qinj, dt.qinj, sig.qinj] = gen_data_in('sihi_smooth(sigadd_quad(\i_inj_x,\i_inj_y))',HitTree);
        int_itor = interp1(t.itorsm, sig.itorsm, tbase);
        int_iquad = interp1(t.qinj, sig.qinj, tbase);
        i_ratio = abs(int_itor./int_iquad);
        %Get average:
        s_ratio = timeBounds(1)*1e-3;
        e_ratio = timeBounds(2)*1e-3;
        [aa, Is] = min((tbase - s_ratio).^2);
        [aa, Ie] = min((tbase - e_ratio).^2);
        CurrAmp = mean(i_ratio(Is:Ie))
        
        if Gain
            x = CurrAmp;
        else
            x = freq(i)-(.5-.25*j);
        end
        
        % get displacement ( vel * interval )
        interval = (1/(2*freq(i)))*100;
        Displacement = (dat.vel -Offset(i)).*interval;
        
        % get average value % quadrature error
        meanDisp = mean(mean(Displacement));
        stdDisp = mean(std(Displacement)); % removed Fitting error: 14.5 doesnt have it yet.
        Error = sqrt( stdDisp.^2 );%+ (mean(mean(dat.velU*interval))).^2 +(mean(std(dat.velU*interval))).^2);
        
        %plot the mean with errorbar
        handle(i) = errorbar(ax,x,meanDisp,Error,Marker{i},'color',ColorVector(j,:),'LineWidth',3);
    end
end

%14.5 Comparison section
timeBounds = [TBase(3).time(1) TBase(3).time(2)];
% Get Current amplification
clear t dt sig;
HitTree.openTree('hitsi',129817);% open shot
% Get data
[t.itorsm, dt.itorsm,sig.itorsm] = gen_data_in('sihi_smooth(\i_tor_spaavg)',HitTree);
[t.qinj, dt.qinj, sig.qinj] = gen_data_in('sihi_smooth(sigadd_quad(\i_inj_x,\i_inj_y))',HitTree);
int_itor = interp1(t.itorsm, sig.itorsm, tbase);
int_iquad = interp1(t.qinj, sig.qinj, tbase);
i_ratio = abs(int_itor./int_iquad);
%Get average:
s_ratio = timeBounds(1)*1e-3;
e_ratio = timeBounds(2)*1e-3;
[aa, Is] = min((tbase - s_ratio).^2);
[aa, Ie] = min((tbase - e_ratio).^2);
CurrAmp = mean(i_ratio(Is:Ie))

for i = 1: size(shots(3).shots,2)

    if Gain
        x = CurrAmp;
    else
        x = freq(3)-(.5-.25*i);
    end
    
    interval = (1/(2*freq(3)))*100;
    Displacement = (shots(3).shots(:,i) -Offset(3)).*interval;
    
    
    meanDisp = mean(Displacement);
    stdDisp = mean(Std(:,i)*interval);
    
    handle(3) = errorbar(ax,x,meanDisp,stdDisp,Marker{3},'color',ColorVector(i,:),'LineWidth',3);
end

if Gain
    XLim = [2 3.5];
    XLab = 'Current Gain ';
    Unit = 'I_{torr} / I_{quad}';
else
    XLim = [0 75];
    XLab = 'Frequency ';
    Unit = '[kHz]';
end

if Pos
    Curr = ', Positive Current';
else
    Curr = ', Negative Current';
end
set(ax,'XLim',XLim);
xlabel([XLab Unit]);
ylabel('Displacement  [cm]');
title(['Average displacement Vs ' XLab Curr]);
%annotation('textbox',[.806,.703,.065,.16],'String',{'14.5 kHz: \diamond';'53.5 kHz: \star';'68.5 kHz: \sq'});
%set(ax,'LineWidth',2);
if Pos
    legend([handle(3) handle(1)],'14.5 kHz','53.5 kHz')
else
    legend([handle(3) handle(1) handle(2)],'14.5 kHz','53.5 kHz','68.5 kHz')
end

%% Displacement scatter plot, as a function of current amplification,
%  Compairing 14.5, 53.5, 63.5 kHz, with error bars
%  - Pull in data, average velocity, multiply by interval.
%  - Get std, get phase, get current amp. Plot.

clear all; close all;
import MDSplus.*

%Load 14.5 velocity separated by phase:
load('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data\Temp\Phase129810.mat');
Velocity = PhaseVelocity.Velocity;
Std = PhaseVelocity.Std;
Std = [Std(:,2) Std(:,4) Std(:,10)];

% 53&68 shots to use ( will find 14)
shots(1).shots = [128585,128594,128592];% THE ORDER IS IMPORTANT:
shots(2).shots = [129214,129213,129215];% IT MAKES PHASE-COLORS ~MATCH
%shots(3).shots = [Velocity(:,2),Velocity(:,4),Velocity(:,10)]; % HARD CODED ORDER: MANUAL SELECTION
shots(3).shots = [Velocity(:,2) Velocity(:,4) Velocity(:,10)];
freq = [53.5, 68.5, 14.5];

% find the zero point, to counteract the offset:
% WARNING: THIS MIGHT NOT BE SCIENCE
% Calculated by averaging the most extreme datapoints.
% For 53, this may be off by 25%
Offset68 = 8.14306718748471; % average of 68(2) and 68(3)
Offset53 = -3.27336872971015;% value for 53(3) ( This makes some sense, based on the relative phases and displacements )
Offset14 = -51.6082215996631; % mean of two colum-means, five columns appart ( approximately the same, regardless of column choice)
Offset = [Offset53 Offset68 Offset14];

TimeBounds = [1.26 1.81]; % curr. gain plateu


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

% Calculate the Current Ampflication
HitTree = Connection('landau.hit');
HitTree.openTree('hitsi',128592);% open shot

ts = -0.0001; % start time for new time base
%ts2 = ts - 0.0002; % start time for cropping data in TEST LATER
te = 0.005; % end time for new time base
%te2 = te + 0.0002; % end time for cropping data in
dtn = 0.4e-6; % dt for new time base (2.5 MHz)
tbase = ts: dtn: te; % new time base

[t.itorsm, dt.itorsm,sig.itorsm] = gen_data_in('sihi_smooth(\i_tor_spaavg)',HitTree);
[t.qinj, dt.qinj, sig.qinj] = gen_data_in('sihi_smooth(sigadd_quad(\i_inj_x,\i_inj_y))',HitTree);
int_itor = interp1(t.itorsm, sig.itorsm, tbase);
int_iquad = interp1(t.qinj, sig.qinj, tbase);
i_ratio = abs(int_itor./int_iquad);
    
s_ratio = TimeBounds(1)*1e-3;
e_ratio = TimeBounds(2)*1e-3;
    
[aa, Is] = min((tbase - s_ratio).^2);
[aa, Ie] = min((tbase - e_ratio).^2);

plot(tbase(Is:Ie).*1e3, i_ratio(Is:Ie)); hold on;
plot(tbase(Is:Ie).*1e3,ones(length(Is:Ie))*mean(i_ratio(Is:Ie)))
figure;
    
% colors
ColorVector=hsv(3);

%h1 = figure('Name','Displacement Comparison'); hold on;
ax = axes('FontSize',18,'LineWidth',3); hold on; grid on;

for i = 1:2; % loop through frequencies
    for j = 1: length(shots(i).shots)
        timeBounds = TimeBounds; % refresh time bounds
        % get the dat from the file
        load(['T:\IDS\Data Repository\TEMP\dat' num2str(shots(i).shots(j)) '10.mat']);
        
        %trim to the channel and time window we want
        while any(any(isnan(dat.vel)))
            dat = trimRandT(dat,chan_range,timeBounds);
            timeBounds(2)=timeBounds(2)-.01;
        end 
        % get displacement ( vel * interval )
        interval = (1/(2*freq(i)))*100;
        Displacement = dat.vel .*interval -Offset(i);
        
        % get average value % quadrature error
        meanDisp = mean(mean(Displacement));
        stdDisp = mean(std(Displacement)); % removed Fitting error: 14.5 doesnt have it yet.
        Error = sqrt( stdDisp.^2 );%+ (mean(mean(dat.velU*interval))).^2 +(mean(std(dat.velU*interval))).^2);
        
        %plot the mean with errorbar
        t = errorbar(ax,freq(i)-(.5-.25*j),meanDisp,Error,'marker','o','color',ColorVector(j,:),'LineWidth',3);
    end
end

%14.5 Comparison section

for i = 1: size(shots(3).shots,2)
    interval = (1/(2*freq(3)))*100;
    Displacement = shots(3).shots(:,i).*interval -Offset(3);
    
    
    meanDisp = mean(Displacement);
    stdDisp = mean(Std(:,i)*interval);
    
    x = errorbar(ax,freq(3)-(.5-.25*i),meanDisp,stdDisp,'marker','o','color',ColorVector(i,:),'LineWidth',3);
end


set(ax,'XLim',[0, 70]);
xlabel('Frequency  [kHz]');
ylabel('Displacement  [cm]');
title('Average displacement Vs Frequency');
%set(ax,'LineWidth',2);
        
%% OKAY. KYLE. FINE.
%clear all;
addpath(genpath('C:\Users\HITSI\Documents\GitHub\IDS_RIAN\NewCodes\'));
in(1).shot = 160609202;
in(1).NIMROD=1;
in(1).line=1;
in(1).timeScale = 1e-3;
in(1).color = {'b',[66, 188, 244]./255};
in(1).style = {'-','--'};
in(1).doubleplot = 1;
in(1).velShift = -5;
in(1).legend={'NIMROD'};
in(2).shot = 160728013;
in(2).NIMROD=0;
in(2).line=2;
in(2).timeScale = 1e-3;
in(2).doubleplot = 1;
in(2).color = {'r',[225,105,0]./255};
in(2).style = {'-','--'};
in(2).legend={'HIT-SI3'};
in(2).velShift = -5;


flipLoImpact=[47].*(in(1).shot>229499);
plotError=0;
lnwdth=2;


%% Parameters


%% Load Data
figure; ax1=axes('Parent',gcf);
for n = 1:length(in)
    if in(n).shot >= 160728009 && in(1).shot <= 160728024
    timebound = [1.65,2.00];
%     chan_ranget = [8:24];
%     chan_rangep = [44:60];
    chan_ranget = [10:32];
    chan_rangep = [37:59];
elseif in(n).shot >= 160609200 in(1).shot <= 160610000 % NIMROD Data
    timebound = [1.6,1.8];
   chan_ranget = [10:32];
    chan_rangep = [37:59];
    xlim=[-20,60];
    
end
chan_range = [chan_ranget, chan_rangep];

if in(n).NIMROD ==0
    input=load(['T:\IDS\Data Repository\dat' num2str(in(n).shot) '10.mat']); % Real HIT-SI Data
else
    input=load(['T:\IDS\Data Repository\NIMROD\dat' num2str(in(n).shot) '10.mat']); % Real HIT-SI Data
end
dat=input.dat;

% Flip Data
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
%error('ABORT\n')
% Trim Data
display(['Pre-Trim Data Length: ' num2str(length(dat(1).impacts))]);
dat = trimRange(dat, chan_range, plotError,timebound.*(1./in(n).timeScale),[]); % for some reason, this wont save to workspace
display(['Post-Trim Data Length: ' num2str(length(dat(1).impacts))]);
assignin('base','datTrim',dat);

%try;saveDat(n).VelError = dat(in(n).line).velU;end
%try;saveDat(n).TempError = dat(in(n).line).tempU;end


%% Data Analysis
dat(in(n).line).vel = averageNans(dat(in(n).line).vel)+in(n).velShift; % remove nans
dat(in(n).line).temp = averageNans(dat(in(n).line).temp); % remove nans
                
if ~isempty(in(n).doubleplot)
    % Find where each fiber bundle begins and ends.
    doubleplot(1,:) = 1:(length(dat(1).impacts))/2;
    doubleplot(2,:) = (length(dat(1).impacts)/2)+1:length(dat(1).impacts);

    % initialize Sine_Fit parameters
    param(:,1,n) = dat(1).impacts(doubleplot(1,:));
    param(:,6,n) = dat(1).impacts(doubleplot(2,:));
end



% Loop Through Channels
for i = 1:length(doubleplot)
    dataAvg(i,1) = mean(dat(in(n).line).vel(:,doubleplot(1,i)));
    dataAvg(i,2) = mean(dat(in(n).line).vel(:,doubleplot(2,i)));
end

% Old Way
data(1:length(dat(1).time),:) = dat(in(n).line).vel(:,1:(length(dat(1).impacts))/2)+in(n).velShift;
data(length(dat(1).time)+1:2*length(dat(1).time),:) = ...
  dat(in(n).line).vel(:,(length(dat(1).impacts)/2)+1:end)+in(n).velShift;
%assignin('base','data',data);
 % find index corresponding to time bounds
nTimeLim(1) = find(dat(1).time.*in(n).timeScale >= timebound(1), 1);
nTimeLim(2) = find(dat(1).time.*in(n).timeScale <= timebound(end), 1, 'last');

for i = 1:2
for m = 1:length(doubleplot);         
    selection = data((nTimeLim(1):nTimeLim(2))+(i-1)*size(data,1)/2, m);
    dataAvg1(m,i) = mean(selection(~isnan(selection)));
    dataStd1(m,i) = std(selection(~isnan(selection)));
    dataDispl1(m,i) = dataStd1(m,i).*(1/14500)./(2*pi) .*1e5;
end
end
%% Plotting
plot(ax1,dat(1).impacts(1:length(dat(1).impacts)/2),-(dataAvg(:,1)-dataAvg(:,2))./2,'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1},...
                             'MarkerEdgeColor',[in(n).color{1}]); hold on;
plot(ax1,dat(1).impacts(1:length(dat(1).impacts)/2),-(dataAvg1(:,1)-dataAvg1(:,2))./2,'color',[in(n).color{2}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{2},...
                             'MarkerEdgeColor',[in(n).color{2}]); hold on;
end

legend(ax1,{ [in(1).legend{1} ': New']  , [in(1).legend{1} ': Old'] , [in(2).legend{1} ': New'] ,[in(2).legend{1} ': Old'] });
title('Toroidal Flow Profile');
xlabel('Impact Parameters [cm]');
ylabel('Velocity [km/s]');
grid on;
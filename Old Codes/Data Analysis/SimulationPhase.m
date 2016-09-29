% Finding velocity by phase for simulation data, where the data points
% arent necessarilly evenly spaced

% find first two times where the sine wave is at zero, linspace 11 points
% between them ( ten and then one again ). Find the two times closest on
% each side to these points. polyfit between each channel, and pull the
% time that we want. 
close all;
clear all;
%sineTimeBounds = [.7068,.93];% time limits for looking at sine graph
sineTimeBounds = [.6635,.8849];
cd('T:\IDS\Data Repository');
addpath('T:\IDS\Display');
shot = [9129499];
cycles = 3; % complete injector cycles

%chan = [50:58]; % spheromak zone
chan = [8:24]; % torroidal region
debug = 1;
load(['dat' num2str(shot) '10.mat']);

[~,minPoint] = min((dat.iinjxTime-sineTimeBounds(1)).^2);
[~,maxPoint] = min((dat.iinjxTime-sineTimeBounds(2)).^2);

dat.iinjxTime = dat.iinjxTime(minPoint:maxPoint);
dat.iinjx = dat.iinjx(minPoint:maxPoint);

dat = trimRange(dat,chan); %trim to the correct range
dat.vel(28,9)=7;

%find the times of the zeros
zeroCross = zeros(cycles+1,1); % times where sine wave starts, needs to see the end of the last cycle
toggle = 0; % sine wave crosses zero twice
itr=1; % zerocross iterator
if debug
    figure;plot(dat.iinjxTime,dat.iinjx); hold on;
    title('X\_Injector with Zeros and Phases');
    xlabel('Time, [ms]');
    ylabel('KV');
end

for i = 1:length(dat.iinjx)-1
    
    % check if zero is between i and i+1, if so, interpolate to find zero
    if ((dat.iinjx(i)>0 && dat.iinjx(i+1)<0) || (dat.iinjx(i)<0 && dat.iinjx(i+1)>0)) && toggle ==0
        [p,s] = polyfit([dat.iinjxTime(i),dat.iinjxTime(i+1)],[dat.iinjx(i),dat.iinjx(i+1)],1);
        y = polyval(p,dat.iinjxTime(i):.0001:dat.iinjxTime(i+1));
        [~,zeroPoint] = min(y.^2); % find where zero happens
        zeroCross(itr) = dat.iinjxTime(i) +.0001*(zeroPoint-1); % find when it happens
        itr=itr+1;
        toggle = 1;
        if debug
            plot(zeroCross(itr-1),y(zeroPoint),'*','color','red');
        end
    elseif ((dat.iinjx(i)>0 && dat.iinjx(i+1)<0) || (dat.iinjx(i)<0 && dat.iinjx(i+1)>0))
        toggle = 0;
    end
end
    

% find the times of the other phases
PhaseTimes = zeros(cycles, 10);
colors = hsv(10);
for i = 1:cycles % zerocross is of size cycles+1, so this is ok.
    Times = linspace(zeroCross(i),zeroCross(i+1),11); 
    PhaseTimes(i,1:10)=Times(1:10);
    if debug
        for j = 1:10
            plot([PhaseTimes(i,j),PhaseTimes(i,j)],[-20,20],'color',colors(j,:));
        end
    end
end


% cycle through the times, find the closes velocity times, linspace all
% over the place, pull the vel points at the time closest to the sine time.
% Store these in the phase array, these should be the correct cords.
Velocity = zeros(size(dat.vel,2),10,cycles);
for i = 1:cycles
    for j = 1:10 % this should get the indexes of the closest time points.
        [Y,I] = sort( (dat.time-PhaseTimes(i,j)).^2 );
        velocityTemp = zeros(size(dat.vel,2), 1);
        for k = 1:size(dat.vel,2) % loop through channels at that time point, find the velocities
            p = polyfit([dat.time(I(1)),dat.time(I(2))],[dat.vel(I(1),k),dat.vel(I(2),k)],1);
            velocityTemp(k) = polyval(p,PhaseTimes(i,j));
        end
        Velocity(:,j,i) = velocityTemp(:);
    end
end


if debug
    colors = hsv(cycles+1);
    for i = 1:10
        figure; hold on;
        for j = 1:cycles
            plot(dat.impacts,squeeze(Velocity(:,i,j)),'color',colors(j,:))
        end
        title(['Phase number: ' num2str(i)]);
        xlabel('Impacts');
        ylabel('Km/s');
        plot(dat.impacts,median(squeeze(Velocity(:,i,:)),2),'-o','color',colors(cycles+1,:));
    end
end


%output
PhaseVelocity.Velocity = median(Velocity,3);
Phase = linspace(.00001,2,11)'; % cant be zero, or comparephase eleminates it
PhaseVelocity.Phase = Phase(1:10);
PhaseVelocity.Impacts = dat.impacts;
PhaseVelocity.Std = std(Velocity,0,3);% technically this is wrong: std is for mean, not median.
cd('T:\IDS\Analysis Repository\Phase Data');
save(['Phase' num2str(shot) ], 'PhaseVelocity');

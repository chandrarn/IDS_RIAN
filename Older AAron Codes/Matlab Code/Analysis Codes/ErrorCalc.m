% ALLL THE ERRROR
% Rian Chandra
% ( requires lm method in batch correct 2 )
% This should add in quadrature all the errors to get a line which
% represents the +- ammount

% load in shot ( in a loop )
% use all timepoints between min and max bounds, to get less error/
% get vel error
% line up the start and end times that compairphase uses EXACTLY
% average across the channels for that time limit
% get std for each channel. 
% add std with average velocity in quad
% Take mean of velocity, get STD, add in quad too
% BOOM

close all; clear all;
torPlot = 1; % true for toroidal array (fibers 1 to 36)
chan_ranget = [8:24];
chan_rangep = [47:61]; % poloida
if torPlot
    chan_range = chan_ranget;
else
    chan_range = chan_rangep;
end
        
shots = [128585, 128592, 128594, 128580, 128581, 128586, 128587, 128588, 128593, 128595, 128596];

%shots = [129213, 129214, 129215];
%timeBounds = [1.0 1.97; 1.0 1.97; 1.0 1.9];
timeBounds = [1.0 1.97;1.0 1.5; 1.0 1.97; 1.0 1.97;1.0 1.97; 1.0 1.97; 1.0 1.97; 1.0 1.97; 1.0 1.97; 1.0 1.97; 1.0 1.97];

% shots = [129810];
% timeBounds = [1.49 1.87;];
shift=0;
saveShot=0;
quad=zeros(16,10);

% load('Phase129810');
% Phase=PhaseVelocity.Phase;
for i = 1:length(shots);
    cd('T:\IDS\Data Repository\TEMP');
    eval(sprintf('load(''dat%i10'');',shots(i)));
    cd('T:\IDS\Display');
    dat = trimRange(dat, chan_range);
    cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data\Temp');
    minPoint = find(dat.time>timeBounds(1));
    minPoint = minPoint(1)+shift+i;
    maxPoint = find(dat.time>timeBounds(2));
    maxPoint = maxPoint(1);
    while find(isnan(dat.velU(minPoint:maxPoint,:)))
        maxPoint=maxPoint-1; % assume no NaN at start
        if maxPoint == 0
            break;
        end
    end
    minPoint
    maxPoint
%     while mod(maxPoint-minPoint,10)~=0 % FOR 14.5 DATA ONLY
%         maxPoint= maxPoint-1;
%     end

    points = minPoint:maxPoint; %for 53, 68
%     points = minPoint:10:maxPoint;
    
    meanVelU = mean(dat.velU(points,:))';
    stdVel = std(dat.vel(points,:))';
    stdVelU = std(dat.velU(points,:))';

    
    quad(:,i)= sqrt( (meanVelU.^2) + (stdVelU.^2) + (stdVel.^2) );
%     figure;
%     plot(quad(:,i));title(['Quadrature Error ' num2str(shots)]);xlabel('Channel Number');ylabel('Km/s');
     figure;
    plot(mean(dat.velU(points,:))','LineWidth',8,'color','green');
    hold on; 
    plot(std(dat.velU(points,:))','LineWidth',8,'color','red');
    plot(dat.velU(points,:)');
    title('Fitting Error');xlabel('Channel Number');ylabel('Km/s');legend('Mean','Std');
    hold off;
    figure; 
    plot(mean(dat.vel(points,:)),'LineWidth',8,'color','green');
    hold on; 
    plot(std(dat.vel(points,:)),'LineWidth',8,'color','red');
    plot(dat.vel(points,:)');
    title('Summation Error');xlabel('Channel Number');ylabel('Km/s'); legend('Mean','Std');
    %cd('T:\IDS\Data Repository\TEMP')
    %temp=quad(:,i);
    hold off

       figure
       plot(quad(:,i),'color','yellow','LineWidth',7); 
       hold on; 
       i
       plot( 1:16,ones(16,1).*mean(quad(:,i),1), 'LineWidth', 7,'color','black');
       
     xlabel('Channel Number');ylabel('Km/s');title(['Phase: ' num2str(shots(i))]);
     hold off
    if saveShot
        temp=quad(:,i);
        save(['quad' num2str(shots) 'Phase ' num2str(shots(i)) '.mat' ], 'temp');
        
%         figure; plot(quad(:,i)); hold on; plot( mean(quad(:,i),1), 'LineWidth', 7);
%         xlabel('Channel Number');ylabel('Km/s');title(['Phase: ' num2str(Phase(i))]);
%         save('quad14','temp');   
    end
     
            
end
cd('T:\IDS\Data Repository\TEMP');
% save('quad14.mat','quad');

    
    
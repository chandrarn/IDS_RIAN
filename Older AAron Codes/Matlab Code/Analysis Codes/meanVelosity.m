% script to plot the mean velocity from dat.Vel, with error bars, to get an
% idea of what the correct offset for the data should be.
% Supplementary code to Compare_Plots 3
close all;

addpath('T:\IDS\Data Repository'); % where corrected IDS data is kept
cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data\Temp')
shot = 12982010;

%eval(sprintf('load(''dat%i'');', shot)); % Real HIT-SI Data

chan_range = [8:27];
timeLim = [1.3,1.7];

% Trim all Data for channel range
dat = trimRange(dat, chan_range);
figure;
plot(dat.vel(:,19));
figure; plot(sum(dat.vel(:,:),2)./length(dat.vel));
data=dat.vel(108:208,:);
%129810 210 250   150 300
% 129817          150 250
% 129819          150 290
% 129820          120 270  

%figure;plot(data(:,8))

for(i=1:19)
    p = polyfit(1:length(data(:,i)),(data(:,i))',1);
    data(:,i)=data(:,i)-(p(1).*(1:length(data(:,i)))+p(2))';
end
%Get mean for each channel
DatStd=std(data)';
size(data(:,19))
size(1:110)
figure; plot(1:101,data(:,19),1:101,mean(data(:,19)).*ones(101)');
% hold on; plot(linspace(30,-28,22),mean(data(:,1)));
% hold off;
DatMean=mean(data(:,:))';
figure;
errorbar(linspace(30,-28,19),DatMean,DatStd);
xlabel('Impact Parameters BACKWARDS');
ylabel('Velocity, Km/S');
title1=['Shot ' num2str(shot) ' average velocity per channel with Std'];
title([title1]);
% %p = polyfit(1:length(data(:,i)),sum(data(:,i),2)',1);
% figure;plot(data(:,8));
% figure; plot(data(:,1));
% figure; plot(data(:,3));
% figure; plot(data(:,5));
% figure; plot(data(:,7));
% figure; plot(data(:,9));
% figure; plot(data(:,11));
% 
% figure;plot(data(:,8));
% figure; plot(data(:,1));
% figure; plot(data(:,3));
% figure; plot(data(:,5));
% figure; plot(data(:,7));
% figure; plot(data(:,9));
% figure; plot(data(:,11));
% data=cumtrapz(data);
% figure; plot(data(:,1));
% figure; plot(data(:,3));
% figure; plot(data(:,5));
% figure; plot(data(:,7));
% figure; plot(data(:,9));
% figure; plot(data(:,11));
% figure; plot(data(:,13));
% figure; plot(data(:,15));
% figure; plot(data(:,17));
% figure; plot(data(:,19));
% figure; plot(data(:,21));
% 
% 
% 
% 
% 
% 

%% Velocity vs impacts at one specific Phase
% we wont know the phase, but every tenth velocity timepoint is at the same
% one. By starting at the same point, we know we are dealing with the same
% dataset.


%data=dat.vel(146:175,:);% new dataset

for(i=0:9)
    for(j=0:2)
        velocity=data(1+i+j*10,:);
    end
    STD=std(velocity);
    MEAN=mean(velocity,1);
    figure;
    errorbar(MEAN,STD);
    xlabel('Impact Parameters BACKWARDS');
    ylabel('Velocity, Km/S, with 1 Std');
    title1=['Phase: ' num2str(i*.2) ' \Pi radians'];
    title([title1]);
end







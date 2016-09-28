%% script to manually find the peaks fclear all;
% load the shot
shot=[129499];
timeStart=.9;
colorVector=hsv(200);
shift=-20;
close all;
cd('T:\IDS\Data Repository\TEMP');%my data
%cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data')%NSTX data
eval(sprintf('load(''dat%i10'');', shot));
cd('T:\IDS\Display');%my data
dat = trimRange(dat, [8:27]);
dat.velU=dat.velU+shift;%or fitting in batch correct

startpoint=198;
endpoint=274; % from compairplots 3 with timeinms=0;
figure;
hold on;

DATMEAN=mean(real(dat.velU(:,startpoint:endpoint)),2);
plot(linspace(-23,28,19),DATMEAN,'LineWidth',6,'color','red');

for i = startpoint:endpoint
    plot(linspace(-23,28,19),dat.velU(:,i));
end

set(gca,'XLim',[-23,28]);
legend('Mean',['Times: ' num2str(dat.time(startpoint)) 'mS-' num2str(dat.time(endpoint)) 'mS']);
xlabel('Impacts');
ylabel('Km/S');
title(['Velocity for shot: ' num2str(shot)]);

DatSTD=std(real(dat.velU(:,startpoint:endpoint))');
figure;
errorbar(linspace(-23,28,19),DATMEAN',DatSTD);
title1=['Shot ' num2str(shot) ' average velocity per channel with Std'];
title([title1]);
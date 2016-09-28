%% Program to detect systemic instrument error by checking if the post-injector shuttoff velocity is nonzero
% go from timestart till NaNs, take all available data
% Rian Chandra


clear all;
% load the shot
shot=[129215];
timeStart=1.5;
colorVector=hsv(200);


cd('T:\IDS\Data Repository');%my data
%cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data')%NSTX data
eval(sprintf('load(''dat%i10'');', shot));
dat = trimRange(dat, [8:27]);
dat.vel=dat.vel';

startPoint=find(dat.time>=timeStart-.001);
startPoint=startPoint(1)
i=0;
TEMP=isnan(dat.vel(:,startPoint+i))
TEMP=max(isnan(dat.vel(:,startPoint+i)))
%Final=zeros(19,22);
figure;
hold on;
% initial = dat.vel(:,startpoint-40ish +i)
% pull double window stuff from compare phase
% switch between plotting axis in single loop?

S=get(0,'ScreenSize');
h1 = figure('Visible','on','Name','Positive/Negative Comparison','Position',...
    [5 35 S(3)-12 S(4)-110],'Color',[1 1 1]);

%Left Graph
ax(1)=axes('Parent',h1,'Position',[.05,.35,.4,.6],'FontSize',15);
set(gca,'XLim',sort([-23,28]),'LineWidth',2,'YLim',[-10,10]);

ylabel('Velocity, km/s');
xlabel('Impact Parameter');
grid on;
hold on
Final=zeros(19,200);
while (max(isnan(dat.vel(:,startPoint+i)))~=1)&&i<=44
    
    Final(:,i+1)=dat.vel(:,startPoint+i)-20;
    plot(linspace(-23,28,19),Final(end:-1:1,i+1),'Color',colorVector(i+1,:));
    i=i+1;
end
title(ax(1),strcat('Post-Shuttoff',num2str(shot),' Time Points: ',num2str(dat.time(startPoint)),' - ',num2str(dat.time(startPoint+i))),'FontSize',18);

%Right Graph
ax(4)=axes('Parent',h1,'Position',[.55,.35,.4,.6],'FontSize',15);
set(gca,'XLim',sort([-23,28]),'LineWidth',2,'YLim',[-10,10]);
title(ax(4),strcat('Pre-Shutoff',num2str(shot),' Time Points: ',num2str(dat.time(startPoint-10)),' - ',num2str(dat.time(startPoint))),'FontSize',18);
ylabel('Velocity, km/s');
xlabel('Impact Parameter');
grid on;
hold on;

for i = 1:30
    
    Initial(:,i)=dat.vel(:,startPoint+i-30)-20;
    plot(linspace(-23,28,19),Initial(end:-1:1,i),'Color',colorVector(i,:));
    i=i+1;
end

xlabel('Impacts');
ylabel('Km/S');

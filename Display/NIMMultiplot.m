% Multiplot
% expansion for plotting nimrod values and averages alongside 
% Note: need to manually save the figure as .eps, otherwise it saves in
% black and white. You know, because fuck all.

close all;
clear all;

% shot number
Ishot = 8150303018;
Nshot = 150318028;

prefix = [0,8];
saving = 0;
plotCurrents = 0;

 plotType = 1; % Velocity
%  plotType = 2; % Temperature
% plotType = 3; % Intensity


% Set these: ( add :2: if you only want every other line
timebound = [1:104];
%chan_range = 50:58;
chan_range = [1:28];
Nshift = 1.07;  % NIMROD time shifting


%set up figure
S=get(0,'ScreenSize');
analysisHeight = S(4)-110;
figureWidth = (S(3)-12)/3;
colors = [ 0 0 1; 12/255 117/255 0; 1 0 0 ];
h = figure('Visible','on','Name','MULTIPLOT-Lines','Position',...
    [5 35 figureWidth analysisHeight],'Color',[1 1 1]);
ax = axes('Parent',h,'Position',[.075,.08,.8,.85]);


hold on;
%THIS WILL BREAK ON LINUS:
cd('T:\IDS\Data Repository');
addpath('T:\IDS\Display');

%load IDS DAta
cd('T:\IDS\Data Repository');
Idat = importdata(['dat' num2str(Ishot) '10.mat'])
% make surf plot to check time bounds
%figure;surf(Idat.vel); shading interp; view( [ 0 90]);
Idat.time = Idat.time' .*1e-6;
 
% NIMROD DATA

cd('T:\IDS\Data Repository\');
Ndat = importdata(['dat8' num2str(Nshot) '10.mat'])
Ndat.time = Ndat.time + Nshift;
Ndat.ItorTime = Ndat.ItorTime +Nshift;

%remove NaNs
for passes =  1:3 % go through three times to make sure
    for m = 12:28 % chans
        for n = 2:(length(Ndat.time)-1) % timepts
            if isnan(Ndat.vel(n,m))
                assignin('base','n',n);
                assignin('base','m',m);
                Ndat.vel(n,m)=mean([Ndat.vel(n-1,m),Ndat.vel(n+1,m)]);
            end
            if isnan(Ndat.temp(n,m))
                assignin('base','n',n);
                assignin('base','m',m);
                Ndat.temp(n,m)=mean([Ndat.temp(n-1,m),Ndat.temp(n+1,m)]);
            end
            if isnan(Ndat.int(n,m))
                assignin('base','n',n);
                assignin('base','m',m);
                Ndat.int(n,m)=mean([Ndat.int(n-1,m),Ndat.int(n+1,m)]);
            end
        end
    end
end

%figure; plot(Idat.ItorTime,Idat.Itor,'b',Ndat.ItorTime,Ndat.Itor,'g',Idat.time,ones(1,length(Idat.time))*15,'b*',Ndat.time,ones(1,length(Ndat.time))*16,'g*');set(gca,'XLim',[0,5]);

% average onto IDS timebase
NIDS.vel = zeros(length(timebound)-1,size(Ndat.vel,2));
NIDS.temp = zeros(length(timebound)-1,size(Ndat.vel,2));
NIDS.int = zeros(length(timebound)-1,size(Ndat.vel,2));
x = 1; % counter
oldX = 1;
% for i = timebound(2:end)
%     while(x<length(Ndat.time)&&Ndat.time(x)<Idat.time(i))
%         x= x+1;% increment counter till we hit the end of one exposure period
%     end
%     NIDS.vel(i-timebound(1),:) = mean(Ndat.vel(oldX:x,:));
%     NIDS.temp(i-timebound(1),:) = mean(Ndat.temp(oldX:x,:));
%     NIDS.int(i-timebound(1),:) = mean(Ndat.int(oldX:x,:));
%     %figure;hold on; plot(Ndat.vel(oldX:x,:)');plot(NIDS.vel(i-timebound(1),:),'r','linewidth',3); plot(Idat.vel(i-1,:)'+2,'linewidth',3);
% 
%     oldX = x;
% end

NIDS = Ndat;


%dat = trimRange(dat, chan_range);

switch plotType
    case 1
        Ndata = NIDS.vel;
        Idata = Idat.vel+2;% add offset if neccessary
        titles = 'Velocity';

        offset = 20; % 129499 velocity
        sidebar = [ num2str(offset) ' km/s per division'];
        
        %assignin('base','TEMPDATA1',data);
    case 2
        Ndata = NIDS.temp;
        Idata = Idat.temp;
        titles = 'Temperature';

        offset = 50; %129499 temp
        sidebar = [ num2str(offset) ' eV per division'];
    case 3
        Ndata = NIDS.int*1e-2;
        Idata = Idat.int*1e-3;
        titles = 'Intensity';
        sidebar = [ 'arb'];
        offset = 50; %129499 int
end

%Trim data to chan_range
Idata = Idata(:,int32(linspace(1,168,28)));
Ndata = Ndata(:,chan_range);

% if(i==1)
%     data = data; % put the IDS offset here if it needs one
%     timeOffset = 0;
% elseif(i==2)
%     timeOffset =1.335; % sometimes Psi-tet doesnt need a time offset
% elseif(i==3)
%     timeOffset = .9377;
% end

% dat.time = dat.time + timeOffset;

% offset each line
for j = 1:size(Idata,2)
   %PhaseVelocity.Velocity(j,:) = PhaseVelocity.Velocity(j,:)+ (j-1)*50;   
    Ndata(:,j) = Ndata(:,j)+ (j-1)*offset;
    Idata(:,j) = Idata(:,j)+ (j-1)*offset;
end


%t(i,:) = plot(ax,0:.2:1.8,PhaseVelocity.Velocity','color',colors(i,:));
%size(dat.time);

t(1,:) = plot(ax,Idat.time(timebound),Idata(timebound,:),'color',colors(1,:),'LineWidth',1.5);
t(2,:) = plot(ax,Idat.time(timebound(1:end)),Ndata(timebound,:),'color',colors(2,:),'LineWidth',1.5);
    
    
for i = 1:size(Idata,2)
    % draw the text
    %y = PhaseVelocity.Velocity(i,10) + 5;
    y = offset*(i-1) +5;
    %text(1.85,y,num2str(PhaseVelocity.Impacts(i)));
    text(Idat.time(timebound(end))+.005,y,num2str(Idat.impacts(i),2),'fontsize',14);
    plot(ax,[Idat.time(1) Idat.time(end)],zeros(2)+(i-1)*offset,'-','color','black');
    %plot(ax,[timebound(1)+.0005 timebound(end)],zeros(2)+(i-1)*offset,'*','color','black');
    hold on;
    
end

% cosmetics
% Legend. THIS IS ACTUALLY KINDA COOL
 %legend(ax,[t(1,1),t(2,1),t(3,1)],'Location','NorthEastOutside',{'IDS','PSI-TET','NIMROD'},'fontsize',14);
legend(ax,[t(1,1),t(2,1)],'Location','BestOutside',{'IDS','NIMROD'},'fontsize',14);
% legend(ax,t(2,1),'Location','BestOutside',{'PSI-TET'});
% legend(ax,t(3,1),'Location','BestOutside',{'NIMROD'});

xlabel('Time [ms]','fontsize',14);
ylabel(sidebar,'fontsize',14);

title(titles,'fontsize',14);
set(ax,'XLim',[ Idat.time(timebound(1)) Idat.time(timebound(end))]);
set(gca,'LineWidth',1.5);
set(gca,'fontsize',14);
box on;
grid on;
set(gca,'YLim',[-25 offset*size(Idata,2)-offset/2]);%-offset/2
set(gca,'YTick',[]);

i = text(Idat.time(timebound(end))+.02,(offset*size(Idata,2))/2,'R  [cm]','fontsize',14);
set(i,'rotation',-90)


if plotCurrents
    % plot currents
    j = figure('Visible','on','Name','MULTIPLOT-Currents','Position',...
        [5 35 figureWidth analysisHeight],'Color',[1 1 1]);
    ax2 = axes('Parent',j,'Position',[.15 .08 .8 .15],'FontSize',15); % currents
    %load(['dat' num2str(Ishots) '10.mat']);

    plot(ax2,Ndat.ItorTime,Ndat.Itor,'LineWidth',1.5); hold on;
    plot(ax2,Idat.ItorTime,Idat.Itor,'color','red','LineWidth',1.5);

    set(gca,'LineWidth',1.5);
    xlabel('Time [ms]');
    ylabel('I_{INJ} [kA]');
    set(ax2,'XLim',[Idat.time(timebound(1)-2) Idat.time(timebound(end)+1)]);
    grid on;
end

%saving
if saving
    cd('T:\IDS\Analysis Repository');
    saveas(h,['Multiplot ' titles num2str(shot(1))],'eps');
    saveas(j,['MultiplotCurr ' num2str(shot(1))],'eps');
end
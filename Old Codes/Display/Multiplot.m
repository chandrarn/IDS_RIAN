% Multiplot
% Note: need to manually save the figure as .eps, otherwise it saves in
% black and white. You know, because fuck all.

%close all;
clear all;

% shot number
shot = 129810;

prefix = [0,5,8];
saving = 0;
plotCurrents = 0;

 plotType = 1; % Velocity
% plotType = 2; % Temperature
% plotType = 3; % Intensity


% Set these: ( add :2: if you only want every other line
timebound = [1.664,1.85];
chan_range = 50:58;
%chan_range = [8:24];

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

for i = 1:3

    if i ==1
        cd('T:\IDS\Data Repository');
        load(['dat' num2str(prefix(i)*1000000 + shot) '10.mat'])

         dat.time = dat.time;% *1e3;
    elseif i == 2 % special psitet velocity case
        cd('T:\IDS\Data Repository\');
        load(['dat' num2str(prefix(i)*1000000 + shot) '10.mat'])

        dat.time = dat.time(1:end);%-0.8298;
    elseif i == 3 % special nimrod temperature case
        cd('T:\IDS\Data Repository\');
        load(['dat' num2str(prefix(i)*1000000 + shot) '10.mat'])


    end

    dat = trimRange(dat, chan_range);
    
    switch plotType
        case 1
            data = dat.vel;
            titles = 'Velocity';
            
            offset = 40; % 129499 velocity
            sidebar = [ num2str(offset) ' km/s per division'];
            if i ==1
                data = data-12;
                assignin('base','TEMPDATA',data);
            end
            assignin('base','TEMPDATA1',data);
        case 2
            data = dat.temp;
            titles = 'Temperature';
            
            offset = 100; %129499 temp
            sidebar = [ num2str(offset) ' eV per division'];
        case 3
            data = dat.int;
            titles = 'Intensity';
            sidebar = [ 'arb'];
            offset = 900; %129499 int
            if i == 1 % if ids, add additional offset
                data = data.*2;
            end
    end
    
    
    if(i==1)
        data = data; % put the IDS offset here if it needs one
        timeOffset = 0;
    elseif(i==2)
        timeOffset =1.335; % sometimes Psi-tet doesnt need a time offset
    elseif(i==3)
        timeOffset = .9377;
    end
    
    dat.time = dat.time + timeOffset;
    
    % offset each line
    for j = 1:size(data,2)
       %PhaseVelocity.Velocity(j,:) = PhaseVelocity.Velocity(j,:)+ (j-1)*50;   
        data(:,j) = data(:,j)+ (j-1)*offset;
    end
    
    
    %t(i,:) = plot(ax,0:.2:1.8,PhaseVelocity.Velocity','color',colors(i,:));
    size(dat.time);
    size(data');
    t(i,:) = plot(ax,dat.time,data,'color',colors(i,:),'LineWidth',1.5);
    
end
    
for i = 1:size(data,2)
    % draw the text
    %y = PhaseVelocity.Velocity(i,10) + 5;
    y = offset*(i-1) +5;
    %text(1.85,y,num2str(PhaseVelocity.Impacts(i)));
    text(timebound(end)+.005,y,num2str(dat.impacts(i),2),'fontsize',14);
    plot(ax,[dat.time(1) dat.time(end)],zeros(2)+(i-1)*offset,'-','color','black');
    %plot(ax,[timebound(1)+.0005 timebound(end)],zeros(2)+(i-1)*offset,'*','color','black');
    hold on;
    
end

% cosmetics
% Legend. THIS IS ACTUALLY KINDA COOL
 legend(ax,[t(1,1),t(2,1),t(3,1)],'Location','NorthEastOutside',{'IDS','PSI-TET','NIMROD'},'fontsize',14);
%legend(ax,[t(1,1),t(2,1)],'Location','NorthEastOutside',{'IDS','PSI-TET'},'fontsize',14);
% legend(ax,t(2,1),'Location','BestOutside',{'PSI-TET'});
% legend(ax,t(3,1),'Location','BestOutside',{'NIMROD'});

xlabel('Time [ms]','fontsize',14);
ylabel(sidebar,'fontsize',14);

title(titles,'fontsize',14);
set(ax,'XLim',timebound);
set(gca,'LineWidth',1.5);
set(gca,'fontsize',14);
box on;
grid on;
set(gca,'YLim',[-25 offset*size(data,2)-offset/2]);%-offset/2
set(gca,'YTick',[]);

i = text(timebound(end)+.02,(offset*size(data,2))/2,'R  [cm]','fontsize',14);
set(i,'rotation',-90)


if plotCurrents
    % plot currents
    j = figure('Visible','on','Name','MULTIPLOT-Currents','Position',...
        [5 35 figureWidth analysisHeight],'Color',[1 1 1]);
    ax2 = axes('Parent',j,'Position',[.15 .08 .8 .15],'FontSize',15); % currents
    load(['dat' num2str(shots(1)) '10.mat']);

    plot(ax2,dat.iinjxTime,dat.iinjx,'LineWidth',1.5); hold on;
    plot(ax2,dat.iinjyTime,dat.iinjy,'color','red','LineWidth',1.5);

    set(gca,'LineWidth',1.5);
    xlabel('Time [ms]');
    ylabel('I_{INJ} [kA]');
    set(ax2,'XLim',timebound);
    grid on;
end

%saving
if saving
    cd('T:\IDS\Analysis Repository');
    saveas(h,['Multiplot ' titles num2str(shot(1))],'eps');
    saveas(j,['MultiplotCurr ' num2str(shot(1))],'eps');
end

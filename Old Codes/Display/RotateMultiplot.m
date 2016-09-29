% Multiplot
% expansion for plotting nimrod values and averages alongside 
% Note: need to manually save the figure as .eps, otherwise it saves in
% black and white. You know, because fuck all.

close all;
clear all;
clc;

% shot number
Ishot =150401027;
Nshot = 150303018;

prefix = [0,8];
saving = 0;
plotCurrents = 0;
plotAverages = 1;

 plotType = 1; % Velocity
%  plotType = 2; % Temperature
% plotType = 3; % Intensity


% Set these: ( add :2: if you only want every other line
timebound = [7:11];
%chanrange = 50:58;
chanrange = [12:27];
Nshift = 1.07;%1.07;  % NIMROD time shifting


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
% figure;
% surf(Idat.vel); shading interp; view( [ 0 90]);

Idat.time = Idat.time' .*1e-6;
 
% NIMROD DATA

cd('T:\IDS\Data Repository\');
Ndat = importdata(['dat8' num2str(Nshot) '10.mat'])
Ndat.time = Ndat.time + Nshift;
%Ndat.ItorTime(1167) = mean([Ndat.ItorTime(1166),Ndat.ItorTime(1168)]);
Ndat.ItorTime = Ndat.ItorTime +Nshift;
Ndat.impacts = Ndat.impacts(end:-1:1);



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
for i = timebound(2:end)
    while(x<length(Ndat.time)&&Ndat.time(x)<Idat.time(i))
        x= x+1;% increment counter till we hit the end of one exposure period
    end
    
    NIDS.vel(i-timebound(1),:) = mean(Ndat.vel(oldX:x,:));
    NIDS.temp(i-timebound(1),:) = mean(Ndat.temp(oldX:x,:));
    NIDS.int(i-timebound(1),:) = mean(Ndat.int(oldX:x,:));
    % MASSIVE CHEATING ABOUT TO OCCUR
    if i ~= 3
        NIDS.vel(i-timebound(1),24) = mean([NIDS.vel(i-timebound(1),23),NIDS.vel(i-timebound(1),25)]);
        NIDS.temp(i-timebound(1),24) = mean([NIDS.temp(i-timebound(1),23),NIDS.temp(i-timebound(1),25)]);
        NIDS.int(i-timebound(1),24) = mean([NIDS.int(i-timebound(1),23),NIDS.int(i-timebound(1),25)]);
    end
    if plotAverages ~=0
        figure;
        %hold on; plot(Ndat.impacts,Ndat.temp(oldX:x,:)');h1=plot(Ndat.impacts,NIDS.temp(i-timebound(1),:),'r','linewidth',3); h2=plot(Idat.impacts(chanrange),Idat.temp(i-1,chanrange)'+2,'linewidth',3);
        %hold on; plot(Ndat.impacts,Ndat.int(oldX:x,:)');h1=plot(Ndat.impacts,NIDS.int(i-timebound(1),:),'r','linewidth',3); h2=plot(Idat.impacts(chanrange),Idat.int(i-1,chanrange)'+2,'linewidth',3);
        hold on; plot(Ndat.impacts-14,Ndat.vel(oldX:x,:)'+7);h1=plot(Ndat.impacts-14,NIDS.vel(i-timebound(1),:)+7,'color',[ 0 .5 0],'linewidth',3); h2=plot(Idat.impacts(chanrange),Idat.vel(i-1,chanrange)'+2,'-*','linewidth',3);

        title(['Time : ' num2str(Idat.time(i-1)) ' [ms]'],'fontweight','bold','fontsize',13);
        grid on;
        legend([h1 h2],{'NIMROD','IDS'});
        set(gca,'Fontweight','Bold');
        set(gca,'fontsize',13);
        box on;
        xlabel('Impact parameter [cm]');
        ylabel('Velocity [km]');
    end
    oldX = x;
end
cd('T:\RChandra\MGH2015\');
if plotAverages == 0


    %dat = trimRange(dat, chanrange);

    switch plotType
        case 1
            Ndata = NIDS.vel;
            Idata = Idat.vel+2;% add offset if neccessary
            titles = 'Velocity';

            offset = 10; % 129499 velocity
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

    %Trim data to chanrange AND TIME RANGE ( chan range may be unnessary now)

    Idata = Idata(timebound(1:(end-1)),chanrange); assignin('base','TIDAT',Idata);
    Ndata = Ndata(:,:);

    % if(i==1)
    %     data = data; % put the IDS offset here if it needs one
    %     timeOffset = 0;
    % elseif(i==2)
    %     timeOffset =1.335; % sometimes Psi-tet doesnt need a time offset
    % elseif(i==3)
    %     timeOffset = .9377;
    % end

    % dat.time = dat.time + timeOffset;

    % offset each line IN TIME, NOT IN SPACE
    for j = 1:size(Ndata,1) % how many time points
       %PhaseVelocity.Velocity(j,:) = PhaseVelocity.Velocity(j,:)+ (j-1)*50;   
        Ndata(j,:) = Ndata(j,:)+ (j-1)*offset;
        Idata(j,:) = Idata(j,:)+ (j-1)*offset;
    end

    assignin('base','TI2',Idata);
    %t(i,:) = plot(ax,0:.2:1.8,PhaseVelocity.Velocity','color',colors(i,:));
    %size(dat.time);
    %size(data');
    t(1,:) = plot(ax,Idat.impacts(chanrange),Idata,'-*','color',colors(1,:),'LineWidth',1.5);
    t(2,:) = plot(ax,Ndat.impacts(:)-14,Ndata+7,'-*','color',colors(2,:),'LineWidth',1.5);


    for i = 1:size(Idata,1)
        % draw the text
        %y = PhaseVelocity.Velocity(i,10) + 5;
        y = offset*(i-1) -5;
        %text(1.85,y,num2str(PhaseVelocity.Impacts(i)));
        %text(Idat.impacts(chanrange(1))+1,y,num2str(Idat.time(i+timebound(1)-1),4),'fontsize',14);
        text(26,y,num2str(Idat.time(i+timebound(1)-1),4),'fontsize',14);
        plot(ax,[Idat.impacts(chanrange(1)) Idat.impacts(chanrange(end))],zeros(2)+(i-1)*offset,'-','color','black');
        %plot(ax,[timebound(1)+.0005 timebound(end)],zeros(2)+(i-1)*offset,'*','color','black');
        hold on;

    end

    % cosmetics
    % Legend. THIS IS ACTUALLY KINDA COOL
     %legend(ax,[t(1,1),t(2,1),t(3,1)],'Location','NorthEastOutside',{'IDS','PSI-TET','NIMROD'},'fontsize',14);
    legend(ax,[t(1,1),t(2,1)],'Location','BestOutside',{'IDS','NIMROD'},'fontsize',14);
    % legend(ax,t(2,1),'Location','BestOutside',{'PSI-TET'});
    % legend(ax,t(3,1),'Location','BestOutside',{'NIMROD'});

    xlabel(ax,'R [cm]','fontsize',14);
    ylabel(ax,sidebar,'fontsize',14);

    title(ax,titles,'fontsize',14);
    %set(ax,'XLim',[ Idat.impacts(chanrange(end)) Idat.impacts(chanrange(1))]);
    set(ax,'XLim',[-25 25]);
    set(ax,'LineWidth',1.5);
    set(ax,'fontsize',14);
    set(ax,'box', 'on');
    grid on;
    set(ax,'YLim',[-7 offset*size(Idata,1)-offset/2]);%-offset/2
    set(ax,'YTick',[]);
    %i = text(Idat.impacts(chanrange(1))+7,(offset*size(Idata,1))/2,'Time  [ms]','fontsize',14);
     i = text(30,(offset*size(Idata,1))/2 +2,'Time  [ms]','fontsize',14);
    set(i,'rotation',-90)


    if plotCurrents
        % plot currents
        j = figure('Visible','on','Name','MULTIPLOT-Currents','Position',...
            [5 35 figureWidth 300],'Color',[1 1 1]);
        ax2 = axes('Parent',j,'Position',[.12 .2 .8 .7],'FontSize',15); % currents
        %load(['dat' num2str(Ishots) '10.mat']);

        p1=plot(ax2,Ndat.ItorTime,Ndat.Itor,'color',[0 .5 0],'LineWidth',1.5); hold on;
        p2=plot(ax2,Idat.ItorTime,Idat.Itor,'color','blue' ,'LineWidth',1.5);

        plot(ax2,Idat.time,ones(1,length(Idat.time))*10,'*','color','blue');
        plot(ax2,Ndat.time,ones(1,length(Ndat.time))*13,'*','color', [0 .5 0]);
        set(gca,'LineWidth',1.5);
        xlabel('Time [ms]');
        ylabel('I_{INJ} [kA]');
        set(ax2,'XLim',[Idat.time(timebound(1)-2) Idat.time(timebound(end)+1)]);
        legend([p1,p2],'NIMROD','IDS');
        grid on;
    end

    %Plot all data for visualization purposes


    %saving
    if saving
        cd('T:\IDS\Analysis Repository');
        saveas(h,['Multiplot ' titles num2str(shot(1))],'eps');
        saveas(j,['MultiplotCurr ' num2str(shot(1))],'eps');
    end
end
cd('T:\RChandra\MGH2015\');
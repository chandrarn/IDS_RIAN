% RIAN CHANDRA MAY 2014
% plot two windows, for negative and positive , with 
%Load in the shots


%close all;
clear all;

% Note: if compairing phases betweem 14 and 53, put 53 on left
%shots = [128580:128581, 128585:128588,128592:128596];
%shots = [129810,129817];
shotsLeft = [129810]; % Positive Current ( verified )
shotsRight = [129817]; % Negative current
%shotsLeft= [128585,128592,128594]; %positive current
%shotsRight= [128580:128581, 128586:128588, 128593,128595:128596];
%shotsLeft = [129810];
%shotsRight = [129819];
%shotsLeft = [129213:129215]; % 68 kHz data Negative

%shotsRight = [];

%settings
velShiftLeft = 0;%-20; % velocity shift ammount
velShiftRight =0;%-20;
chansLeft=19;% number of channels on left
chansRight=19;
titleLeft = 'Positive Toroidal Current';
titleRight= 'Negative Toroidal Current';
currentLeft = [];
currentRight = [];
Is14=1; % 14.5kHz shot, changes display system.
IsCompare=0; % 14.5 shot on right to compair to 53 on left.
ColorVector=hsv(10); % max number of colors to display in the vector
PhaseWing = .2;
saveFigure = 1; % save figure to file
fileName = ['T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data\Temp\' num2str(shotsLeft) 'VS' num2str(shotsRight)]; % file name for .png image
PlotRight = 1; % If only plotting one graph
Neg53Vec=[1,0,0.600000000000000;1,0,0;0.800000000000001,0,1;0,0.400000000000000,1;1,0.600000000000000,0;0,1,1;;0.199999999999999,0,1;0,1,0.400000000000000;1,0,0;];
Pos53Vec=[1,0.600000000000000,0;1,0,0.600000000000000;0.800000000000000,1,0;];
Neg68Vec=[0.200000000000000,1,0;1,0,0;0.800000000000001,0,1;];

%cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data');
cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data\Temp');
LeftCount=1;
RightCount=1;

VelLeft=zeros(chansLeft,10);%Note: cannot display more than 10 phases at once;
VelRight=zeros(chansRight,10);
PhaLeft=zeros(10,10);
PhaRight=PhaLeft;
%Notes for compairing phase
%Where phasevelocity is 14.5 and phase1 is from 53
% find phasevelocity.phase > phase1 -phasewing && phasevelocity.phase
% <phase1 + phasewing
% check to see which is close: if multiple: wont be more than two as long
% as phase wing is less than .2 radians
% probably set it to like .15. Technically anything .1 and greater will
% always work ( directly in between phases )
% Put velocity and phase from PhaseVelocity at that point in new list,
% correspoinding to chronology of 53
% can use shotsRigt and left for loading. Need to load both shots
% completely, parts of 14.5 will be eliminated in a for loop, for size 10,
% afterwards.
% check if phase(i) is withing phasewing, if so, add to list, if not, move
% on. 


for(i=shotsLeft)
    eval(sprintf('load(''Phase%i'');', i));
    PhaseVelocity.Velocity=PhaseVelocity.Velocity+velShiftLeft;
    if(Is14)
             PhaLeft(:,1)=PhaseVelocity.Phase;%if phase is vector, if 14.5
             %size(Vel
             VelLeft(:,:)=(PhaseVelocity.Velocity);
    else
             PhaLeft(LeftCount,1)=PhaseVelocity.Phase;%in case Phase is scalar
             VelLeft(:,LeftCount)=(PhaseVelocity.Velocity);
    end
    LeftCount=LeftCount+1;
    ImpactsLeft=PhaseVelocity.Impacts;%this will overwrite every time, but
                                 %it should be the same value, so its ok.
end

 
for(i=shotsRight)
    TEMp=i
    clear PhaseVelocity;
    eval(sprintf('load(''Phase%i'');', i));
    PhaseVelocity.Velocity=PhaseVelocity.Velocity+velShiftRight;
    if(Is14||IsCompare)
             PhaRight(:,1)=PhaseVelocity.Phase;%if phase is vector, if 14.5
             VelRight(:,:)=PhaseVelocity.Velocity;
    else

             PhaRight(RightCount,1)=PhaseVelocity.Phase;%in case Phase is scalar
             VelRight(:,RightCount)=PhaseVelocity.Velocity;
    end
    RightCount=RightCount+1;
    ImpactsRight=PhaseVelocity.Impacts;
end

 %eliminate zeroed columns and rows, if 53.5 and less data
    PhaLeft( ~any(PhaLeft,2), : ) = [];  %rows
    PhaLeft( :, ~any(PhaLeft,1) ) = [];  %columns
    PhaRight( ~any(PhaRight,2), : ) = [];  
    PhaRight( :, ~any(PhaRight,1) ) = []; 
    VelLeft( ~any(VelLeft,2), : ) = []; 
    VelLeft( :, ~any(VelLeft,1) ) = [];  
    VelRight( ~any(VelRight,2), : ) = [];  
    VelRight( :, ~any(VelRight,1) ) = [];  
    
sinewave=sin(0:pi/100:2*pi);

%In case there are too many shots to display
if(length(shotsLeft)>5)
    shotsLeft=shotsLeft(1:5);
end
if(length(shotsRight)>5)
    shotsLeft=shotsRight(1:5);
end


%% Compair phase 14.5 vs 53.5
% this should cycle through 14.5 and compair all the phases to a given
% phase from 53, to find the one that matches most closely, and put it in
% the corresponding spot on the list, so colors will be same.
if(IsCompare)
    newPhaRight=zeros(length(PhaLeft),1);
    newVelRight=zeros(length(VelRight),length(PhaLeft));
    for(i=PhaLeft')
        for(j=PhaRight')
            if( (i-PhaseWing) < j && j < (i+PhaseWing) ) % if ja withing phasewing
                newPhaRight(find(PhaLeft==i),1)=j;
                newVelRight(:,find(PhaLeft==i))=VelRight(:,find(PhaRight==j));
            end
        end
    end
    length(newPhaRight)
    PhaRight=zeros(length(PhaLeft),1);
    VelRight=zeros(length(newVelRight),length(PhaRight));
    PhaRight=newPhaRight; % reassign
    VelRight=newVelRight;
end


%% PLOTTING SECTION

%make the window
S=get(0,'ScreenSize');
h1 = figure('Visible','on','Name','Positive/Negative Comparison','Position',...
    [5 35 S(3)-12 S(4)-110],'Color',[1 1 1]);

%First Sine plot
ax(1)=axes('Parent',h1,'Position',[.05 .08 .4 .19],'FontSize',15);
plot(0:1/100:2,sinewave,'b',0:1/100:2,-cos(0:pi/100:2*pi),'r','LineWidth',2);
legend('X-Inj','Y-Inj');
set(gca,'LineWidth',2);
hold on;
%plot(0:1/100:2,-cos(0:pi/100:2*pi),'LineWidth',2,'Color','red');
hold on
for(i=1:length(PhaLeft))
    %[abs(sin( (.4*(i+2)/length(PhaLeft))*2*pi)),abs(sin( (.4*(i)/length(PhaLeft))*2*pi  )),abs(sin( (.1*(i)/length(PhaLeft))*2*pi))]
    plot([PhaLeft(i),PhaLeft(i)],[1,-1],'Color',ColorVector(i,:),'LineWidth',2);
end
xlabel(strcat('Phase [degrees]'),'FontSize',20);
ylabel({'Injector';'Current [Arb.]'},'FontSize',20);
set(gca,'YTickLabel',[]);
set(gca,'XTick',linspace(0,2,9));
set(gca,'XTickLabel',[0:45:360]);

if PlotRight
%Second Sine Plot
ax(2)=axes('Parent',h1,'Position',[.55 .08 .4 .19],'FontSize',15);
plot(0:1/100:2,sinewave,0:1/100:2,-cos(0:pi/100:2*pi),'r','LineWidth',2);
legend('X-Inj','Y-Inj');
set(gca,'LineWidth',2);
hold on
for(i=1:length(PhaRight))
    plot([PhaRight(i),PhaRight(i)],[1,-1],'Color',ColorVector(i,:),'LineWidth',2);
end
xlabel(strcat('Phase [degrees]'),'FontSize',20);
ylabel({'Injector';'Current [Arb.]'},'FontSize',20);
set(gca,'YTickLabel',[]);
set(gca,'XTick',linspace(0,2,9));
set(gca,'XTickLabel',[0:45:360]);
end

%Left Graph
ax(3)=axes('Parent',h1,'Position',[.05,.35,.4,.6],'FontSize',15);
set(gca,'XLim',sort([ImpactsLeft(1),ImpactsLeft(end)]),'LineWidth',2,'YLim',[-10,10]);
title(ax(3),[titleLeft ],'FontSize',20);
ylabel('Velocity [km/s]','FontSize',20);
xlabel('Impact Parameter [cm]','FontSize',20);
%shotsLeft= [128585,128592,128594]; %positive current
%shotsRight= [128580:128581, 128586:128588, 128593,128595:128596];

grid on;
box on;
hold on;

for(i=1:length(PhaLeft))
    %[abs(sin( (.4*(i+2)/length(PhaLeft))*2*pi)),abs(sin( (.4*(i)/length(PhaLeft))*2*pi  )),abs(sin( (.1*(i)/length(PhaLeft))*2*pi))]
    plot([ImpactsLeft],squeeze(VelLeft(:,i)),'Color',ColorVector(i,:),'LineWidth',2);
end
%h=legend('Shot: 129213','Shot: 129214','Shot: 129215');
%get(h,'title')
%set(get(h,'title'),'string','Shot #:');

if PlotRight
%Right Graph
ax(4)=axes('Parent',h1,'Position',[.55,.35,.4,.6],'FontSize',15);
set(gca,'XLim',sort([ImpactsRight(1),ImpactsRight(end)]),'LineWidth',2,'YLim',[-10,10]);
title(ax(4),[titleRight],'FontSize',20);
ylabel('Velocity [km/s]','FontSize',20);
xlabel('Impact Parameter [cm]','FontSize',20);

grid on;
box on;
%set(gca,'FontSize',20);
hold on
for(i=1:length(PhaRight))
    plot(linspace(ImpactsRight(1),ImpactsRight(end),length(VelRight)),squeeze(VelRight(:,i)),'Color',[ColorVector(i,:)],'LineWidth',2);
end
%h=legend('Shot: 128580','Shot: 128581', 'Shot: 128586','Shot: 128587','Shot: 128588', 'Shot: 128593','Shot: 128595','Shot: 128596');
%set(get(h,'title'),'string','Shot #:','FontSize',15);
end

%% Save

if saveFigure
    fig_save = getframe(h1);
    [Xfig, mapfig] = frame2im(fig_save);
    imwrite(Xfig, [fileName '.png']);
end

% RIAN CHANDRA MAY 2014
% Loads two shots, plots them with cooresponding colors signifying
% cooresponding phases. Gets phase data from LoopPhase->VelocityPhase. Can
% perform rudementary quantatative analysis of displacement and velocity
% for 14.5 kHz shots. Can compair between 14.5 kHz and 53/68 at
% approximately the same phase. 
% Added analysis plot: can find a basic linear rotation velocity, and a
% basic bulk modulus motion. Note: does not work well for highly nonlinear
% regions of the machiene, such as looking into the injector.


%NOTES:
% 129499 POLOIDAL SHALL HENCEFORTH BE KNOWN AS 1294990
% 1294990 and 81294990 45 degree port real vs nimrod
% 129518 and 8129518 270 degree port (y-inj) real vs nimrod
% 129530 and 8129530 same
% 1129499 and 8129499 torroidal midplane, real vs nimrod
% 1129499 and 5129499 torroidal midplane, real vs psitet
% 1294990 and 51294990 45 degree port, real vs psitet
% 129530 and 5129530 270 degree port (y-inj), real vs psitet
% 129810 and 5129810 135 degree port, real vs psitet


%close all;
clear all;

% Note: if compairing phases betweem 14 and 53, put 53 on left
%shots = [128580:128581, 128585:128588,128592:128596];
%shots = [129810,129817];
%shotsLeft = [129810];
%shotsRight = [129817]; %Positive current
%shotsRight= [128585,128592,128594]; %positive current
%shotsLeft= [128580:128581, 128586:128588, 128593,128595:128596];
%shotsLeft = [129810];
%shotsRight = [129819];
shotsLeft = [129530];
shotsRight = [];

%% SETTINGS
velShiftLeft =0;%-20; % velocity shift ammount
velShiftRight = 0;%-20;
chansLeft=7% number of channels on left
chansRight=16;
titleLeft = 'IDS Data';
titleRight= '14.5 kHz, PSI-TET Data, shot: ';
currentLeft = [' Positive'];
currentRight = [' Positive'];
Is14=1; % 14.5kHz shot, changes display system.
IsCompare=0; % 14.5 shot on right to compair to 53 on left.
ColorVector=hsv(10); % max number of colors to display in the vector
PhaseWing = .2;
save = 0;
plotAnalysis = 1; % plot averaged data, etc. Only works for 14.5
autoCurrent = 0; % get current sign from tree, only works for 14.5, real data


%cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data');
cd('T:\IDS\Analysis Repository\Phase Data');
LeftCount=1;
RightCount=1;

% Get plotting settings
S=get(0,'ScreenSize');
% check if we're plotting square or not 
if isempty(shotsRight)
    figureWidth = (S(3)-12)/3;
    plotWidth = .8;
    plotLeft = .1;
    analysisHeight = 3*(S(4)-110)/4;
    dataHeight = 3*(S(4)-110)/4;
else
    figureWidth = (S(3)-12);
    plotWidth = .4;
    plotLeft = .05;
    analysisHeight = S(4)-110;
    dataHeight = S(4)-110;
end

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

%% Load Shots
for(i=shotsLeft)
    eval(sprintf('load(''Phase%i'');', i));
    PhaseVelocity.Velocity=PhaseVelocity.Velocity+velShiftLeft;
    if(Is14)
             PhaLeft(:,1)=PhaseVelocity.Phase;%if phase is vector, if 14.5
             size(PhaseVelocity.Velocity) % keep this unsurpressed: its important
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

    

%In case there are too many shots to display
if(length(shotsLeft)>5)
    shotsLeft=shotsLeft(1:5);
end
if(length(shotsRight)>5)
    shotsLeft=shotsRight(1:5);
end

%% Generate Current sign from tree
if autoCurrent
    import MDSplus.*
    HitTree = Connection('landau.hit'); % use connection object. Cause
    sig = '\i_tor_spaavg';
    numb = [shotsLeft, shotsRight];
    name = ['Left', 'Right'];
    for i = 1:2
    %left plot
        HitTree.openTree('hitsi',numb(i));
        Torr = NATIVEvalue(HitTree.get(sig)); % get time and signal, find what the signal is
        dt = NATIVEvalue(HitTree.get(['samplinginterval(' sig ')'])); % at a time in the middle.
        tmin = NATIVEvalue(HitTree.get(['minval(dim_of(' sig '))'])); % Check sign of signal.
        tlength = NATIVEvalue(HitTree.get(['size(dim_of(' sig '))']));
        t = tmin + dt*(double(1: tlength) - 1);
        [Y,I] = min( (t-.0015).^2 ); % the current at 1.5mS shouldn't have collapsed yet
        if Torr(I) <0
            eval(['current' name(i) ' = "Negative"']);
        else
            eval(['current' name(i) ' = "Positive"']);
        end
    end
    
%     %right plot
%     HitTree.openTree('hitsi',shotsRight);
%     Torr = NATIVEvalue(HitTree.get(sig));
%     dt = NATIVEvalue(HitTree.get(['samplinginterval(' sig ')']));
%     tmin = NATIVEvalue(HitTree.get(['minval(dim_of(' sig '))']));
%     tlength = NATIVEvalue(HitTree.get(['size(dim_of(' sig '))']));
%     t = tmin + dt*(double(1: tlength) - 1);
%     [Y,I] = min( (t-.0015).^2 );
%     if Torr(I) <0
%         currentRight = ' Negative';
%     else
%         currentRight = ' Positive';
%     end
end
    

%% Analyze the velocity and displacement data
if plotAnalysis
    h = figure('Visible','on','Name','Positive/Negative Comparison','Position',...
    [5 35 figureWidth analysisHeight],'Color',[1 1 1]);
    % get quantititave values for velocity, and error ( only if plotAnalysis )
    plotL = 0; % plot left first
    for i = {'Left','Right'}
        i = char(i);
        if isempty(shotsRight) && isequal(i,'Right')
            break;
        end
        % initialize
        eval(['Velocity = Vel' i ';'])
        eval(['Impacts = Impacts' i ';'])
        eval(['Phase = Pha' i ';'])
        eval(['Slope' i '=zeros(10,1);'])
        eval(['SlErr' i '=zeros(10,1);'])
        eval(['Displ' i ' = zeros(size(Velocity,1),1);'])
        eval(['Diffr' i ' = zeros(size(Velocity,1),1);'])
        for j = 1:10 % get velocity values and error
            [P,F] = polyfit(Impacts,Velocity(:,j),1);
            Y = polyval(P,Impacts,F);
            eval(['Slope' i '(j) = (Y(end)-Y(1))/2;'])
            eval(['SlErr' i '(j) = std(Velocity(:,j)-Y);'])
        end
        % get displacement values and average change in displacement
        for j = 1:size(Velocity,1) 
            eval(['Displ' i '(j) = std(Velocity(j,:));'])
            eval(['Diffr' i '(j)= mean(abs(diff(Velocity(j,:))));'])    
        end

        %Plotting section:
        %Velocity
        ax(2*plotL+1)=axes('Parent',h,'Position',[(plotLeft+.5*plotL),.6,plotWidth,.35],'FontSize',15);
        %set(gca,'XLim',[eval(['Pha' i '(1)']),eval(['Pha' i '(end)'])],'LineWidth',2);% removed YLim
        set(gca,'XLim',[0 2],'LineWidth',2);
        set(gca,'YLim',[-25 20]);
        title(ax(2*plotL+1),['Rotation'],'FontSize',18); % use to have side title
        ylabel('Velocity [km/s]');
        xlabel('Phase');
        %set(ax(2*plotL+1),'YTick',linspace(-5,5,5));
        grid on;
        hold on;
        errorbar(eval(['Pha' i]),eval(['Slope' i]),eval(['SlErr' i]),'lineWidth',2);
        %Displacement
        ax(2*plotL+2)=axes('Parent',h,'Position',[(plotLeft+.5*plotL),.10,plotWidth,.35],'FontSize',15);
        set(gca,'XLim',[eval(['Impacts' i '(1)'])-1,eval(['Impacts' i '(end)'])+1],'LineWidth',2);
        set(gca,'YLim',[0 15]);
        title(ax(2*plotL+2),['Displacement'],'FontSize',18);
        ylabel('Displacement [mm]'); % km/s * 1e-3S exposure * 1000mm/m = meters
        xlabel('Impacts');
        grid on;
        hold on;
        errorbar(eval(['Impacts' i]),eval(['Displ' i]),eval(['Diffr' i]),'lineWidth',2);
        plotL=1;
    end
    %save
    cd('T:\IDS\Analysis Repository');
    saveas(h,['Analysis ' num2str(shotsLeft) 'AND' num2str(shotsRight)],'bmp');
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
h1 = figure('Visible','on','Name','Positive/Negative Comparison','Position',...
    [5 35 figureWidth dataHeight],'Color',[1 1 1]);

sinewave=sin(0:pi/100:2*pi); % generate sine wave

%First Sine plot
ax(1)=axes('Parent',h1,'Position',[plotLeft .08 plotWidth .19],'FontSize',15);
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
xlabel(strcat('Current:',currentLeft,'. Phase in pi radians'));

%Second Sine Plot
if ~isempty(shotsRight) % in case we're plotting square
ax(2)=axes('Parent',h1,'Position',[.55 .08 .4 .19],'FontSize',15);
plot(0:1/100:2,sinewave,0:1/100:2,-cos(0:pi/100:2*pi),'r','LineWidth',2);
legend('X-Inj','Y-Inj');
set(gca,'LineWidth',2);
hold on
for(i=1:length(PhaRight))
    plot([PhaRight(i),PhaRight(i)],[1,-1],'Color',ColorVector(i,:),'LineWidth',2);
end
xlabel(strcat('Current:',currentRight,'. Phase in pi radians'));
end

%Left Graph
ax(3)=axes('Parent',h1,'Position',[plotLeft,.35,plotWidth,.6],'FontSize',15);
set(gca,'XLim',sort([ImpactsLeft(1),ImpactsLeft(end)]),'LineWidth',2,'YLim',[-25,15]);
title(ax(3),titleLeft,'FontSize',18);
ylabel('Velocity, km/s');
xlabel('Impact Parameter');
grid on;
hold on
for(i=1:length(PhaLeft))
    %[abs(sin( (.4*(i+2)/length(PhaLeft))*2*pi)),abs(sin( (.4*(i)/length(PhaLeft))*2*pi  )),abs(sin( (.1*(i)/length(PhaLeft))*2*pi))]
    plot([ImpactsLeft],squeeze(VelLeft(:,i)),'Color',ColorVector(i,:),'LineWidth',2);
end

%Right Graph
if ~isempty(shotsRight) 
ax(4)=axes('Parent',h1,'Position',[.55,.35,.4,.6],'FontSize',15);
set(gca,'XLim',sort([ImpactsRight(1),ImpactsRight(end)]),'LineWidth',2,'YLim',[-18,18]);
title(ax(4),strcat(titleRight,num2str(shotsRight)),'FontSize',18);
ylabel('Velocity, km/s');
xlabel('Impact Parameter');
grid on;
hold on
for(i=1:length(PhaRight))
    % x used to be: linspace(ImpactsRight(1),ImpactsRight(end),length(VelRight))
    plot([ImpactsRight],squeeze(VelRight(:,i)),'Color',[ColorVector(i,:)],'LineWidth',2);
end
end

%saving
if save
    cd('T:\IDS\Analysis Repository');
    saveas(h1,[ num2str(shotsLeft) ' AND ' num2str(shotsRight) ],'bmp');
end


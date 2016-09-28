%plot nimrod magnetics, kinetics, and current
addpath('T:\RChandra');
%load('post_process_output.mat');
close all;
sihi = 1;
xLim =[];% [.55,.925];

%get data
ItorTime = nimsave.time([1:1166 1168:end])*1e3; % that one nimrod datapoint
Itor = nimsave.Bprob.surf.Itor_spa_avg([1:1166 1168:end])*1e-4;
Ia = nimsave.inj.Ia([1:1166 1168:end])*1e-4;
Ib = nimsave.inj.Ib([1:1166 1168:end])*1e-4;
Ic = nimsave.inj.Ic([1:1166 1168:end])*1e-4;
ModeTime = nimsave.energy.et*1e3;
MagModes = nimsave.energy.emn;
KinModes = nimsave.energy.ekn;

%sihi_smooth
ItorS = sihi_smooth(Itor,ItorTime,47.5);
for i = 1:22
    MagModesS(i,:) = sihi_smooth(MagModes(i,:),ModeTime,47.5);
    KinModesS(i,:) = sihi_smooth(KinModes(i,:),ModeTime,47.5);
end

% variables
topOfScreen = .05;
bottomOfScreen = .05;
spacing = .03;
usableArea = 1-topOfScreen-bottomOfScreen-2*spacing;
L1 = bottomOfScreen;
H1 = usableArea*(2/12);
L2 = L1+H1+spacing;
H2 = usableArea*(5/12);
L3 = L2+H2+spacing;
H3 = usableArea*(5/12);
fntsz = 12;
labels = {'I_{Torr} [kA]','Kinetic Modes [J]','Magnetic Modes [J]'};
Ylims = [0,2.5;0,1.25;0,70];


S = get(0, 'ScreenSize');
analysisHeight = S(4) - 110;
figureWidth = (S(3) - 12)/5;
colors = [ 0 0 1; 12/255 117/255 0; 1 0 0 ];
h = figure('Visible', 'on', 'Name', 'NIMROD-Modes', 'Position',...
    [5, 35, figureWidth, analysisHeight], 'Color', [1 1 1]);
ax(1) = axes('Parent', h, 'Position', [0.075, L1, .8, H1]);
ax(2) = axes('Parent', h, 'Position', [0.075, L2, .8, H2]);
ax(3) = axes('Parent', h, 'Position', [0.075, L3, .8, H3]);

%plotting
plot(ax(1),ItorTime,Itor,ItorTime,Ia,ItorTime,Ib,ItorTime,Ic);
plot(ax(2),ModeTime,KinModes);
plot(ax(3),ModeTime,MagModes);
if sihi
    for i =1:3;hold(ax(i),'on');end
    plot(ax(1),ItorTime,ItorS,'k','linewidth',2);
    plot(ax(2),ModeTime,KinModesS,'k','linewidth',2);
    plot(ax(3),ModeTime,MagModesS,'k','linewidth',2);
end

%formatting
linkaxes(ax,'x');
if ~isempty(xLim);set(ax(1),'xlim',xLim);end
for i = 1:3
    set(ax(i),'ylim',Ylims(i,:));
    set(ax(i),'fontweight','bold');
    set(ax(i),'linewidth',2);
    if(i~=1)
        set(ax(i),'xticklabel',[]);
    end
    yl = get(ax(i),'ylim');
    xl = get(ax(1),'xlim');
    set(h,'currentaxes',ax(i));
    t = text(xl(2)+(xl(2)-xl(1))/10,yl(2)*2/3,labels{i},'fontweight','bold');
    set(t,'rotation',-90);
    grid(ax(i),'on');
end
xlabel(ax(1),'Time [mS]','fontweight','bold');
tl =title('NIMROD Modes','fontsize',fntsz,'fontweight','bold');
tlpos = get(tl,'position');
set(tl,'Position',[tlpos(1),tlpos(2)+1,tlpos(3)]);
legend(ax(3),{'\textbf{\textit{n} = 0}','\textbf{\textit{n} = 1}','\textbf{\textit{n} = 2}','\textbf{\textit{n} $\geq$ 3}'},'Location','Best','interpreter','latex');


%figure; plot(ModeTime,MagModes); hold on; plot(ModeTime,MagModesS,'linewidth',2);

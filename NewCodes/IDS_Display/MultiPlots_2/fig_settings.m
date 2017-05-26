%% Figure setup for Multiplot_2
% set up all figures, return handles to axes.
% 

function [h,ax] = fig_settings(plt,in)

h = zeros(12,1);
S = get(0, 'ScreenSize');
analysisHeight = S(4) - 1;
figureWidth = (S(3) - 12)/2.25;
Widen = 150; % 100 or zero usually
colors = [ 0 0 1; 12/255 117/255 0; 1 0 0 ];
h(1) = figure('Visible', 'on', 'Name', ['MULTIPLOT-Lines: ' num2str(in(1).line)], 'Position',...
    [5, 35, figureWidth, analysisHeight], 'Color', [1 1 1]);
colorOrder = get(groot,'defaultAxesColorOrder');
partialColor = {'r',colorOrder(2,:); 'b',colorOrder(1,:);,'g',colorOrder(5,:)};
if ~plt.compactCurrents
    ax(1) = axes('Parent', h(1), 'Position', [0.075, 0.08, 0.8, 0.85]);
else
    ax(1) = axes('Parent', h(1), 'Position', [0.075+(.0585)/2, 0.28, 0.8-.0585, 0.65]);
end
hold(ax(1),'on');

if plt.Averages || isempty(in(1).fftPlot)
    h(2) = figure('Visible', 'on', 'Name', ['MULTIPLOT-Toroidal Flow: ' num2str(in(1).line)], 'Position',...
        [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
    ax(2) = axes('Parent', h(2), 'Position', [0.075, 0.15, 0.85, 0.75]);
    h(4) = figure('Visible', 'on', 'Name', ['MULTIPLOT-Displacement: ' num2str(in(1).line)], 'Position',...
        [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
    hold on;
    grid on;
end
if plt.Averages && ~isempty(in(1).fftPlot)
     h(3) = figure('Visible', 'on', 'Name', ['MULTIPLOT-Phase: ' num2str(in(1).line)], 'Position',...
        [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
    ax(3) = axes('Parent', h(3), 'Position', [0.075, 0.15, 0.85, 0.75]); hold on; grid on;
    h(4) = figure('Visible', 'on', 'Name', ['MULTIPLOT-Displacement: ' num2str(in(1).line)], 'Position',...
        [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
    ax(4) = axes('Parent', h(4), 'Position', [0.075, 0.15, 0.85, 0.75]); hold on; grid on;
     h(5) = figure('Visible', 'on', 'Name', ['MULTIPLOT-FlowRanges: ' num2str(in(1).line)], 'Position',...
        [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
    ax(5) = axes('Parent', h(5), 'Position', [0.075, 0.15, 0.85, 0.75]);
    hold on;
    grid on;
elseif ~plt.Averages && ~isempty(in(1).fftPlot) % This one usually Executes
    %h(2)=figure;
    % Main Analysis Figure
    h(6) = figure('Visible', 'on', 'Name', ['MULTIPLOT-Analysis: ' num2str(in(1).shot)], 'Position',...
    [5, 1+25, figureWidth-450+Widen, analysisHeight-25], 'Color', [1 1 1]);
    ax(6) = axes('Parent',h(6),'Position',[0.15,.08,.8,.25]); hold on; grid on;box on; % Phse
    ax(7) = axes('Parent',h(6),'Position',[0.15,.38,.8,.25]); hold on; grid on;box on; % Flow
    ax(8) = axes('Parent',h(6),'Position',[0.15,.68,.8,.25]); hold on; grid on;box on; % Disp
    mBox=uicontrol('Style','Text');% Make the title
    set(mBox,'String',in(1).AnalysisTitle); set(mBox,'fontweight','bold'); set(mBox,'fontsize',13);
    set(mBox,'Position', [ 30+ceil(Widen/2),960,400,20]);set(mBox,'BackgroundColor',[1 1 1]);
    if plt.Type==1
        title(ax(8),'Toroidal Displacement Phase');
        title(ax(7),'Toroidal Flow');
        title(ax(6),'Maximum Displacement');
    elseif plt.Type==2
        title(ax(8),'Temperature Phase');
        title(ax(7),'Temperature Profile');
        title(ax(6),'Temperature Oscillations');
    end
    if plt.includeTemp
        ax(17)=axes('Parent',h(6),'Position',[0.15,.08,.8,.18]); hold on; grid on;box on;
        title(ax(17),'Temperature Profile');
        set(ax(6),'Position',[0.15,.30,.8,.18]); 
        set(ax(7),'Position',[0.15,.52,.8,.18]);
        set(ax(8),'Position',[0.15,.74,.8,.18]);
    end
    
    % Plot FFT Power Spectrum
    if plt.FFT
        h(7) = figure('Visible', 'on', 'Name', ['MULTIPLOT-FFT Spectrum: ' num2str(in(1).line)], 'Position',...
            [5, 35, figureWidth, 0.35 * analysisHeight], 'Color', [1 1 1]);
        ax(9) = axes('Parent', h(7), 'Position', [0.075, 0.15, 0.85, 0.75]); hold on; grid on;box on;
    end
    
    if plt.SanityPhase
        % Show the phases as calculated various methods
        h10=figure('Visible', 'on', 'Name', ['Phase: Sanity Check ' in(1).legend], 'Position',...
        [5, 1, figureWidth-300, analysisHeight], 'Color', [1 1 1]);
        ax(10) = axes('Parent',h10,'Position',[0.15,.08,.8,.25]); hold on; grid on;box on;
        ax(11) = axes('Parent',h10,'Position',[0.15,.38,.8,.25]); hold on; grid on;box on;
        ax(12) = axes('Parent',h10,'Position',[0.15,.68,.8,.25]); hold on; grid on;box on;
        % Reconstruct a few lines as a sanity check    
        h11=figure('Visible', 'on', 'Name', ['Phase: Reconstruction ' in(1).legend], 'Position',...
        [5, 1, figureWidth-300, analysisHeight], 'Color', [1 1 1]);
        ax(13) = axes('Parent',h11,'Position',[0.15,.08,.8,.25]); hold on; grid on;box on;
        ax(14) = axes('Parent',h11,'Position',[0.15,.38,.8,.25]); hold on; grid on;box on;
        ax(15) = axes('Parent',h11,'Position',[0.15,.68,.8,.25]); hold on; grid on;box on;
    end
    
    if plt.ExplainReconst
        % Data Reconstruction Figure
        h(12)=figure('Visible','on','Name','Explainatory Reconstruction');
        ax(16) = axes('Parent',h(12)); hold on; grid on; box on;
    end
end
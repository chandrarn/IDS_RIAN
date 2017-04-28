%% Settings
clear all; close all; clc;
addpath('~/Magnetics/filters');
addpath('~/Magnetics/general');
addpath('~/Magnetics/BD_SP');
addpath('T:\IDS\Analysis Repository');
import MDSplus.*
doFFT = 0;

shots = [129499, 160728013];

n = 1;
shot(n).title = ['HIT-SI, \Delta' '\phi = 90^{\circ}'];

n = 2;
shot(n).title = ['HIT-SI3, \Delta' '\phi = 120^{\circ}'];

for n = 1:length(shots)
    shot(n).xTick = 0:5:50; % [cm]
    shot(n).xLim = [0, 48]; % [cm]
    
    shot(n).dispLim = [0, 8];
    
    shot(n).tempLim = [0, 48];
    
    shot(n).phaseLim = [-260, 180];
    
    shot(n).flowLim = [-3.6, 6.2];
end

%% Set up Figure
fntsz = 16;
color = ['-k'; '-r'; '-g'; '-b'; '-c'; '-m';
     ':k'; ':r'; ':g'; ':b'; ':c'; ':m'];
S = get(0, 'ScreenSize');

figlabs = ['abcd'; 'efgh'; 'ijkl'; 'mnop'];

dx = 0.45;

dy = 0.211;
y(1) = 0.075;
y(2) = y(1) + dy + 0.01;
y(3) = y(1) + 2*dy + 0.02;
y(4) = y(1) + 3*dy + 0.03;

x = [0.07, 0.54, 0.07, 0.54]; % for big 2x2 arrangement
Y = [0, 0, 0, 0];

ax = [];

H.fig(1) = figure('Color', [1 1 1], 'Position', [0.05 * S(3), 0.05 * S(4), 0.6 * S(3), 0.85 * S(4)]); % currents, Fourier, Inj only, Plasma only

for n = 1:length(shots)
    
    load(['saveDat' num2str(shots(n)) '_4']);
    
    if size(saveDat(1).Displacement, 2) == 2 % check for lower fiber data
        lower = 1;
    else
        lower = 0;
    end
    
    figure(H.fig(1));
    
    %% Displacement
    ax(length(ax) + 1) = subplot('position', [x(n), Y(n) + y(4), dx, dy]);
    set(gca, 'FontSize', fntsz);
    hold on; box on; grid on;
    if shots(n)<229499
        H.p(1) = errorbar(saveDat(1).Impacts, saveDat(1).Displacement(:, 1), saveDat(1).DisplacementError(:, 1), '-r');
        H.p(2) = errorbar(saveDat(2).Impacts, saveDat(2).Displacement(:, 1), saveDat(2).DisplacementError(:, 1), '-b');
    else % hitsi3
        H.p(1) = errorbar(saveDat(1).Impacts, saveDat(1).Displacement(:, 1), saveDat(1).DisplacementError(:, 1), '-r');
        H.p(2) = errorbar(saveDat(2).Impacts, saveDat(2).Displacement(:, 1), saveDat(2).DisplacementError(:, 1), '-b');
        H.p(3) = errorbar(saveDat(1).Impacts, saveDat(1).Displacement(:, 2), saveDat(1).DisplacementError(:, 2), '--r');
        H.p(4) = errorbar(saveDat(2).Impacts, saveDat(2).Displacement(:, 2), saveDat(2).DisplacementError(:, 2), '--b');
    end
    legend(H.p(1:2), ['Shot ' num2str(saveDat(1).shot) ' (+)'], ...
        ['Shot ' num2str(saveDat(2).shot) ' (-)'], 'Location', 'West');
    set(gca, 'XLim', shot(n).xLim, 'XTickLabel', [], 'XTick', shot(n).xTick);
%     text(0.55, 0.06, ['shot ' num2str(shots(n))], 'Units', 'normalized', 'FontSize', fntsz-2);
    text(0.32, 0.9, 'Max Displacement', 'Units', 'normalized', 'FontSize', fntsz);
    if n == 1
        ylabel('[cm]');
    else
        set(gca, 'YTickLabel', []);
    end
    set(gca, 'YLim', shot(n).dispLim);
    text(0.02, 0.9, ['(' figlabs(n, 1) ')'], 'Units', 'normalized', 'FontSize', fntsz);
    
    title(shot(n).title);
    
    %% Phase
    
    % adjust phases to put X or A inj at 0.
    off = saveDat(1).injPhase(1); % offset = X or A inj phase
    saveDat(1).injPhase = saveDat(1).injPhase - off;
    saveDat(1).Phase = saveDat(1).Phase - off;
    saveDat(2).Phase = saveDat(2).Phase - off;
    
    ax(length(ax) + 1) = subplot('position', [x(n), Y(n) + y(3), dx, dy]);
    set(gca, 'FontSize', fntsz);
    hold on; box on; grid on;
    if shots(n)<229499
        H.l(1) = plot(shot(n).xLim, [saveDat(1).injPhase(1), saveDat(1).injPhase(1)], '-m','linewidth',1.5);
        H.l(2) = plot(shot(n).xLim, [saveDat(1).injPhase(2), saveDat(1).injPhase(2)], '-c','linewidth',1.5);
        set(H.l(2), 'Color', 0.8 * [0, 1, 1], 'LineWidth', 1.5);
        legend(H.l, 'X inj', 'Y inj', 'Location', 'Best');
        
        H.p(1) = errorbar(saveDat(1).Impacts, saveDat(1).Phase(:, 1), saveDat(1).PhaseError(:, 1), '-r');
        H.p(2) = errorbar(saveDat(2).Impacts, saveDat(2).Phase(:, 1), saveDat(2).PhaseError(:, 1), '-b');
    else % hitsi3
        H.l(1) = plot(shot(n).xLim, [saveDat(1).injPhase(1), saveDat(1).injPhase(1)], '-m','linewidth',1.5);
        H.l(2) = plot(shot(n).xLim, [saveDat(1).injPhase(2), saveDat(1).injPhase(2)], '-c','linewidth',1.5);
        set(H.l(2), 'Color', 0.8 * [0, 1, 1], 'LineWidth', 1.5);
        H.l(3:4) = plot(shot(n).xLim, [saveDat(1).injPhase(3), saveDat(1).injPhase(3)], '-g', ...
            shot(n).xLim, 360 + [saveDat(1).injPhase(3), saveDat(1).injPhase(3)], '-g','linewidth',1.5);
        set(H.l(3:4), 'Color', [0, 0.7, 0]);
        legend(H.l, 'A inj', 'B inj', 'C inj', 'Location', 'Best');
        
        H.p(1) = errorbar(saveDat(1).Impacts, saveDat(1).Phase(:, 1), saveDat(1).PhaseError(:, 1), '-r');
        H.p(2) = errorbar(saveDat(2).Impacts, saveDat(2).Phase(:, 1), saveDat(2).PhaseError(:, 1), '-b');
        H.p(3) = errorbar(saveDat(1).Impacts, saveDat(1).Phase(:, 2), saveDat(1).PhaseError(:, 2), '--r');
        H.p(4) = errorbar(saveDat(2).Impacts, saveDat(2).Phase(:, 2), saveDat(2).PhaseError(:, 2), '--b');
    end
    set(gca, 'XLim', shot(n).xLim, 'XTickLabel', [], 'XTick', shot(n).xTick);
    text(0.3, 0.9, 'Phase of Displacement', 'Units', 'normalized', 'FontSize', fntsz);
    if n == 1
        ylabel('[deg]');
    else
        set(gca, 'YTickLabel', []);
    end
    set(gca, 'YLim', shot(n).phaseLim);
    text(0.02, 0.9, ['(' figlabs(n, 2) ')'], 'Units', 'normalized', 'FontSize', fntsz);
    
    %% Flow
    ax(length(ax) + 1) = subplot('position', [x(n), Y(n) + y(2), dx, dy]);
    set(gca, 'FontSize', fntsz);
    hold on; box on; grid on;
    plot(shot(n).xLim, [0, 0], '-k');
    if shots(n)<229499
        saveDat(1).FlowError = saveDat(1).FlowError'; % damn it Rian
        saveDat(2).FlowError = saveDat(2).FlowError'; % damn it Rian
        
        H.p(1) = errorbar(saveDat(1).Impacts, saveDat(1).Flow(:, 1), saveDat(1).FlowError(:, 1), '-r');
        H.p(2) = errorbar(saveDat(2).Impacts, saveDat(2).Flow(:, 1), saveDat(2).FlowError(:, 1), '-b');
    else % hitsi3
        saveDat(1).FlowError = saveDat(1).FlowError'; % damn it Rian
        saveDat(2).FlowError = saveDat(2).FlowError'; % damn it Rian
        
        H.p(1) = errorbar(saveDat(1).Impacts, saveDat(1).Flow(:, 1), saveDat(1).FlowError(:, 1), '-r');
        H.p(2) = errorbar(saveDat(2).Impacts, saveDat(2).Flow(:, 1), saveDat(2).FlowError(:, 1), '-b');
    end
    set(gca, 'XLim', shot(n).xLim, 'XTickLabel', [], 'XTick', shot(n).xTick);
    text(0.32, 0.9, 'Net Toroidal Flow', 'Units', 'normalized', 'FontSize', fntsz);
    if n == 1
        ylabel('[km/s]');
    else
        set(gca, 'YTickLabel', []);
    end
    set(gca, 'YLim', shot(n).flowLim);
    text(0.02, 0.9, ['(' figlabs(n, 3) ')'], 'Units', 'normalized', 'FontSize', fntsz);
    
    %% Temperature
    ax(length(ax) + 1) = subplot('position', [x(n), Y(n) + y(1), dx, dy]);
    set(gca, 'FontSize', fntsz);
    hold on; box on; grid on;
    if shots(n)<229499
        saveDat(1).Temp = saveDat(1).Temp'; % damn it Rian
        saveDat(2).Temp = saveDat(2).Temp'; % damn it Rian
        saveDat(1).TempError = saveDat(1).TempError'; % damn it Rian
        saveDat(2).TempError = saveDat(2).TempError'; % damn it Rian
        
        H.p(1) = errorbar(saveDat(1).Impacts, saveDat(1).Temp(:, 1), saveDat(1).TempError(:, 1), '-r');
        H.p(2) = errorbar(saveDat(2).Impacts, saveDat(2).Temp(:, 1), saveDat(2).TempError(:, 1), '-b');
    else % hitsi3
        H.p(1) = errorbar(saveDat(1).Impacts, saveDat(1).Temp(:, 1), saveDat(1).TempError(:, 1), '-r');
        H.p(2) = errorbar(saveDat(2).Impacts, saveDat(2).Temp(:, 1), saveDat(2).TempError(:, 1), '-b');
        H.p(3) = errorbar(saveDat(1).Impacts, saveDat(1).Temp(:, 2), saveDat(1).TempError(:, 2), '--r');
        H.p(4) = errorbar(saveDat(2).Impacts, saveDat(2).Temp(:, 2), saveDat(2).TempError(:, 2), '--b');
    end
    set(gca, 'XLim', shot(n).xLim, 'XTick', shot(n).xTick);
    text(0.36, 0.9, 'Temperature', 'Units', 'normalized', 'FontSize', fntsz);
    xlabel('Chord Impact Parameter [cm]');
    if n == 1
        ylabel('[eV]');
    else
        set(gca, 'YTickLabel', []);
    end
    set(gca, 'YLim', shot(n).tempLim);
    text(0.02, 0.9, ['(' figlabs(n, 4) ')'], 'Units', 'normalized', 'FontSize', fntsz);
    

end

if doFFT
    c = get(0,'DefaultAxesColorOrder');
    
%     figure;
%     ax1=axes('Parent',gcf);
    legCount=1;
    for n = 1:length(shots)
        load(['saveDat' num2str(shots(n)) '_4']);

        if size(saveDat(1).Displacement, 2) == 2 % check for lower fiber data
            lower = 1;
        else
            lower = 0;
        end 
        plot(ax1,saveDat(1).Impacts,saveDat(1).FFT(:,1),'-*','color',c(2*(n-1) +1,:),'linewidth',2);
        hold on;
        leg(legCount)={[num2str(saveDat(1).shot) ': Upper']};legCount=legCount+1;
        plot(ax1,saveDat(2).Impacts,saveDat(2).FFT(:,1),'-*','color',c(2*(n-1) +2,:),'linewidth',2);
        leg(legCount)={[num2str(saveDat(2).shot) ': Upper']};legCount=legCount+1;
        if lower==1
            plot(ax1,saveDat(1).Impacts,saveDat(1).FFT(:,2),'--*','linewidth',2,'color',c(2*(n-1) +1,:));
            leg(legCount)={[num2str(saveDat(1).shot) ': Lower']};legCount=legCount+1;
            plot(ax1,saveDat(2).Impacts,saveDat(2).FFT(:,2),'--*','linewidth',2,'color',c(2*(n-1) +2,:));
            leg(legCount)={[num2str(saveDat(2).shot) ': Lower']};legCount=legCount+1;
        end
    end
    legend(ax1,leg);
    plot(ax1,[saveDat(1).Impacts(1),saveDat(1).Impacts(end)],[.5,.5],'k--');
    grid on; 
    title('Normalized FFT Power Spectrum');
    xlabel('Chord Impact Parameter [cm]');
    ylabel('P_{Inj}/P_{tot}');
    set(gca,'ylim',[0,1]);
    text(.05, 0.9, ['(a)'], 'Units', 'normalized', 'FontSize', fntsz);
end




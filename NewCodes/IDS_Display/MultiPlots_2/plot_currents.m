%% Plot Currents for Multiplots_2
function [] = plot_currents(compactCurrents,ax,h,dat,in,Itor)
    if ~compactCurrents
        j = figure('Visible', 'on', 'Name', 'MULTIPLOT-Currents', 'Position',...
            [5, 35, figureWidth, analysisHeight], 'Color', [1 1 1]);
        ax2 = axes('Parent', j, 'Position', [0.15, 0.08, 0.8, 0.15], 'FontSize', fntsz); % currents
    else
        ax2Pos = [.075+(.0585/2), .08, .8-.0585, .15];
        ax2 = axes('Parent', h, 'Position',ax2Pos , 'FontSize', fntsz);
    end
    if length(in) > 1
        load(['dat' num2str(in(1).shot) '10.mat']);
    end
    if dat(1).shotRef >999999
        plot(ax2, dat(1).iinjaTime.*in(1).injTimeScale, ...
            dat(1).iinja.*in(1).injScale,'LineWidth', lnwdth','color',[0    0.4470    0.7410]); 
        hold on;
        plot(ax2, dat(1).iinjbTime.*in(1).injTimeScale, ...
            dat(1).iinjb.*in(1).injScale, 'color',[0.8500    0.3250    0.0980], 'LineWidth', lnwdth);
        hold on;
        plot(ax2, dat(1).iinjcTime.*in(1).injTimeScale, ...
            dat(1).iinjc.*in(1).injScale, 'color',[0.9290    0.6940    0.1250], 'LineWidth', lnwdth);
        hold on;
    else
         plot(ax2, dat(1).iinjxTime.*in(1).injTimeScale, ...
            dat(1).iinjx.*in(1).injScale,'LineWidth', lnwdth,'color',[0    0.4470    0.7410]); 
        hold on;
        plot(ax2, dat(1).iinjyTime.*in(1).injTimeScale, ...
            dat(1).iinjy.*in(1).injScale, 'color',[0.8500    0.3250    0.0980], 'LineWidth', lnwdth);
        hold on
        
    end
    for i = 1:length(in)
        plot(ax2, dat(1).ItorTime.*in(1).injTimeScale, ...
            Itor(:,i).*in(1).injScale, 'color', in(i).color{1}, 'LineWidth', lnwdth);
    end
    if dat(1).shotRef >999999
        currleg=legend(ax2,{'I_{A}','I_{B}','I_{C}',['I_{torr}'],['I_{torr}']},'location','EastOutside');
    else
        shotName = strsplit(num2str([in(:).shot]));
        currleg=legend(ax2,{'I_{x}','I_{y}',['I_{torr}'],['I_{torr}']},'location','EastOutside');
    end
    set(gca, 'LineWidth', lnwdth);
    set(gca,'fontsize',fntsz);
    set(currleg,'fontsize',fntsz-7);
    set(ax2,'Position',ax2Pos);
    currxl=xlabel('Time [ms]');
    curryl=ylabel('I [kA]');
    currylPos=get(curryl,'Position');
    set(curryl,'Position',[currylPos(1)+.02,currylPos(2:3)]);
    linkaxes([ax,ax2],'x');
    set(ax2, 'XLim', timebound);
    grid on;
end
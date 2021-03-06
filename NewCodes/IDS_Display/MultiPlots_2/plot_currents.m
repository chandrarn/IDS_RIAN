%% Plot Currents for Multiplots_2

function [ax] = plot_currents(compactCurrents,ax,h,dat,in,Itor,ItorTime,fntsz,lnwdth,timebound)

    if ~compactCurrents % Separate currents plot
        j = figure('Visible', 'on', 'Name', 'MULTIPLOT-Currents', 'Position',...
            [5, 35, figureWidth, analysisHeight], 'Color', [1 1 1]);
        ax(2) = axes('Parent', j, 'Position', [0.15, 0.08, 0.8, 0.15], 'FontSize', fntsz); % currents
    else
        ax2Pos = [.075+(.0585/2), .08, .8-.0585, .15];
        ax(2) = axes('Parent', h(1), 'Position',ax2Pos , 'FontSize', fntsz);
    end
    
    % How many injectors?
    if dat(1).shotRef >999999
        %% HIT-SI3
        plot(ax(2), dat(1).iinjaTime.*in(1).injTimeScale, ...
            dat(1).iinja.*in(1).injScale,'LineWidth', lnwdth','color',[0    0.4470    0.7410]); 
        hold on;
        plot(ax(2), dat(1).iinjbTime.*in(1).injTimeScale, ...
            dat(1).iinjb.*in(1).injScale, 'color',[0.8500    0.3250    0.0980], 'LineWidth', lnwdth);
        hold on;
        plot(ax(2), dat(1).iinjcTime.*in(1).injTimeScale, ...
            dat(1).iinjc.*in(1).injScale, 'color',[0.9290    0.6940    0.1250], 'LineWidth', lnwdth);
        hold on;
    else
        %% HIT-SI
         plot(ax(2), dat(1).iinjxTime.*in(1).injTimeScale, ...
            dat(1).iinjx.*in(1).injScale,'LineWidth', lnwdth,'color',[0    0.4470    0.7410]); 
        hold on;
        plot(ax(2), dat(1).iinjyTime.*in(1).injTimeScale, ...
            dat(1).iinjy.*in(1).injScale, 'color',[0.8500    0.3250    0.0980], 'LineWidth', lnwdth);
        hold on
        
    end
    
    % Plot Toroidal current
    for i = 1:length(in)
        plot(ax(2), dat(1).ItorTime.*in(1).injTimeScale, ...
            Itor(:,i).*in(1).injScale, 'color', in(i).color{1}, 'LineWidth', lnwdth);
    end
    
    % Figure and Legend Properties
    if dat(1).shotRef >999999
        currleg=legend(ax(2),{'I_{A}','I_{B}','I_{C}',['I_{torr}'],['I_{torr}']},'location','EastOutside');
    else
        shotName = strsplit(num2str([in(:).shot]));
        currleg=legend(ax(2),{'I_{x}','I_{y}',['I_{torr}'],['I_{torr}']},'location','EastOutside');
    end
    set(gca, 'LineWidth', lnwdth);
    set(gca,'fontsize',fntsz);
    set(currleg,'fontsize',fntsz-7);
    set(ax(2),'Position',ax2Pos);
    currxl=xlabel('Time [ms]');
    curryl=ylabel('I [kA]');
    currylPos=get(curryl,'Position');
    set(curryl,'Position',[currylPos(1)+.02,currylPos(2:3)]);
    linkaxes([ax(1),ax(2)],'x');
    set(ax(2), 'XLim', timebound);
    grid on;
end
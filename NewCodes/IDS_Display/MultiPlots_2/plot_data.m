%% Plot Data for Multiplots_2, with Line Offsets

function [saveDat,t] = plot_data(data,offset,dat,ax,time,errorL,errorU,...
    lnwdth,in,n,saveDat,plotTor,timebound,fntsz,plt,sidebar,titles)

    %% offset each line for plot
    zeroline = zeros(ceil(size(data,1)/2),size(data,2));
    for j = 1:size(data, 2) 
        data(:, j) = data(:, j) + (j-1) * offset;
        zeroline(:,j) = zeros(ceil(size(data,1)/2),1)+(j-1) * offset;
    end

    %% Plot Data
    if in(n).error % if plotting raw data with errobars. Unused.
        time = ndgrid(dat(1).time, 1:size(data, 2));
        t(n, :) = errorbar(ax(1), time, data, errorL, errorU, 'color', in(n).color, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
    else
        if ~isempty(in(n).doubleplot) && ~plotTor % plot both fibers
            saveDat(n).UpperFiberData = data(1:length(dat(1).time),:);
            saveDat(n).LowerFiberData = data(length(dat(1).time)+1:end,:);
            saveDat(n).Impacts = dat(1).impacts(1:size(data,2));
            saveDat(n).Time =time(1:length(dat(1).time));
            if ~in(n).error % Errobars on Data
                t(n, :) = plot(ax(1), time(1:length(dat(1).time)), data(1:length(dat(1).time),:), ...
                    'color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
                t(n, :) = plot(ax(1), time(length(dat(1).time)+1:end)', data(length(dat(1).time)+1:end,:), ...
                    'color', in(n).color{2}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{2});
            else
                 t(n, :) = errorbar(ax(1), time(1:length(dat(1).time)), data(1:length(dat(1).time),:),...
                     dat(in(n).line).velU(:,doubleplot(1,:)),dat(in(n).line).velL(:,doubleplot(1,:)),...
                    'color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
                 t(n, :) = errorbar(ax(1), time(length(dat(1).time)+1:end)', data(length(dat(1).time)+1:end,:), ...
                    dat(in(n).line).velU(:,doubleplot(1,:)),dat(in(n).line).velL(:,doubleplot(1,:)),...
                    'color', in(n).color{2}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{2});
            end 
            plot(ax(1), [time(1),time(length(dat(1).time))], zeroline([1,size(data,1)./2],:), ...
                '--k', 'LineWidth', .5, 'LineStyle', in(n).style{1});
        elseif  ~isempty(in(n).doubleplot) && plotTor % plot the toroidal flow ( UNUSED)
            t(n, :) = plot(ax(1), time(1:length(dat(1).time)), data(length(dat(1).time)+1:end,:)-data(1:length(dat(1).time),:) +zeroline, ...
                'color', 'k', 'LineWidth', lnwdth, 'LineStyle', in(n).style);
            plot(ax(1), [time(1),time(length(dat(1).time))], zeroline([1,size(data,1)./2],:), ...
                '--k', 'LineWidth', .5, 'LineStyle', in(n).style);
        else % HIT-SI, Single Array Plot
             saveDat(n).UpperFiberData = data(1:length(dat(1).time),:);
             saveDat(n).Impacts = dat(1).impacts(1:size(data,2));
            saveDat(n).Time =time(1:length(dat(1).time));
            if ~in(n).error
                t(n, :) = plot(ax(1), time, data, 'color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
            else
                size(data)
                size(dat(in(n).line).velU)
                size(dat(in(n).line).velL)
                t(n, :) = errorbar(ax(1), time, data,dat(in(n).line).velU,dat(in(n).line).velL,...
                'color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
            end
        end
    end
    
    %% Impact Parameter Labels on Right
    for n = 1:size(data, 2)
        y = offset * (n-1) + 0.1 * offset;
        x=timebound(end) + 0.02 * (timebound(end) - timebound(1));
        text(x, y,num2str(dat(1).impacts(n), 2), 'fontsize', fntsz);
        plot(ax(1), [dat(1).time(1), dat(1).time(end)], zeros(2) + (n-1) * offset, '-', 'color', 'k');
    end
    
    %% Misc. Figure Properties
    if ~plt.compactCurrents
        xlabel('Time [ms]','fontsize', fntsz);
    else
        set(gca,'xticklabel',[]);
    end
    
    ylabel(sidebar, 'fontsize', fntsz);
    datatl=title([in(1).AnalysisTitle  ':' titles], 'fontsize', fntsz);
    set(gca, 'LineWidth', lnwdth);
    set(gca, 'fontsize', fntsz);
    box on;
    grid on;
    
    if or(plt.Type == 2, plt.Type == 3) % temperature or Intensity
        yLowerLim = 0;
        yUpperLim = offset * size(data, 2);
    else
        yLowerLim = -0.5 * offset;
        yUpperLim = offset * size(data, 2) - 0.5 * offset;
    end
    set(gca, 'YLim', [yLowerLim, yUpperLim]);
    set(gca, 'YTick', []);

    nt = text(timebound(end) + 0.125 * (timebound(end) - timebound(1)), (offset * size(data, 2))/2, 'Impact Parameter  [cm]', 'fontsize', fntsz);
    set(nt, 'rotation', -90)
end
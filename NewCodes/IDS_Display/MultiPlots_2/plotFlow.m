%% Calculate HIT-SI3 Flow Profiles for Multiplots_2
function [saveDat] = plotFlow(param,dataAvg,flowSupress,plotError,...
    saveDat,ax,h,n,in,plotType,dat,doubleplot,data,RMS,lnwdth)

if isempty(in(1).fftPlot)
    %% If we're not  doing the FFT
    display('NO FFT PLOT')
    L=nanmean(param,2)'; % in this case, param is the STD of the velocity
    U=-nanmean(param,2)';
    figure(h(2)); hold on;
    errorbar(dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2))./2,L,U,'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1});
    drawnow;
    set(gca,'ylim',[-2,20]);

     plot(ax(2),dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2))./2,'color' ,'k','marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1});
    plot(ax(2),xlim,[0,0],'--k')
    ylabel('Toroidal Flow [km/s]'); set(gca,'ylim',[-10,10]);
    xlabel('Impacts [cm]');
    grid on; box on;
    
    saveDat=[];

else
    %% We're doing the FFT
    for i = 1: 1+in(n).doubleplot
        if plotType==1 % Plotting Velocity
             if any(flowSupress(n,:))
                 dataAvg(flowSupress(n,:),1)=NaN;
             end

             % Plotting with Error Bars
             if ~plotError
                 saveDat(n).Flow = -(dataAvg(:,1)-dataAvg(:,2))./2;
                 t3(n)=plot(ax(7),dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2))./2,'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1},...
             'MarkerEdgeColor',[in(n).color{1}]);
             else
                 error = sqrt( mean(dat(in(n).line).velU(:,doubleplot(1,:)).^2) + mean(dat(in(n).line).velU(:,doubleplot(2,:)).^2) + (squeeze(RMS(1,i,n,100))).^2+ (squeeze(RMS(2,i,n,100))).^2 )/sqrt(2*length(dat(1).time));
                 saveDat(n).Flow = -(dataAvg(:,1)-dataAvg(:,2))./2;
                 saveDat(n).FlowError = error;
                 t3(n)=errorbar(ax(7),dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2))./2,error,'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1},...
            'MarkerEdgeColor',[in(n).color{1}]);
             end
        elseif plotType==2 % Plotting Temperature
             t3(n)=plot(ax(7),dat(1).impacts(1:size(data,2)),dataAvg(:,1),'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1},...
            'MarkerEdgeColor',[in(n).color{1}]);
             t3(n)=plot(ax(7),dat(1).impacts(1:size(data,2)),dataAvg(:,2),'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{2},...
            'MarkerEdgeColor',[in(n).color{1}]);
        end

    end
    set(ax(7),'ylim',[-2,20]);
end

end
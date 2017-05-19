%% Calculate HIT-SI3 Flow Profiles for Multiplots_2
function [saveDat] = plotFlow(param,dataAvg,flowSupress,plotError,...
    saveDat,ax,h,n,in,plotType,dat,doubleplot,data,RMS,lnwdth)

for i = 1:size(data, 2)
    cycle1 = (param(i,3,n)*sin(param(i,4,n)+(2*pi).*(0:(1/100):1))+param(i,2,n));
    cycle2 = (param(i,8,n)*sin(param(i,9,n)+(2*pi).*(0:(1/100):1))+param(i,7,n));
    maxDispl(i) = max(cycle1-cycle2);
    minDispl(i) = min(cycle1-cycle2);
end


% If we're not  doing the FFT
if isempty(in(1).fftPlot)
    L=-(dataAvg(:,1)-dataAvg(:,2))+minDispl'; % This may be wrong. We might just want max/min Displ
    U=-(dataAvg(:,1)-dataAvg(:,2))-maxDispl';
    figure(h5)
    errorbar(dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2))./2,L,U,'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style);
    set(gca,'ylim',[-2,20]);

    figure(h3);
    
% We're doing the FFT
else
    if plotType==1 % Plotting Velocity
         if any(flowSupress(n,:))
             dataAvg(flowSupress(n,:),1)=NaN;
         end

         if ~plotError
             saveDat(n).Flow = -(dataAvg(:,1)-dataAvg(:,2))./2;
             t3(n)=plot(ax(7),dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2))./2,'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1},...
         'MarkerEdgeColor',[in(n).color{1}]);
         else
             error = sqrt( mean(dat(in(n).line).velU(:,doubleplot(1,:)).^2) + mean(dat(in(n).line).velU(:,doubleplot(2,:)).^2) + (squeeze(RMS(1,i,n,100))).^2+ (squeeze(RMS(2,i,n,100))).^2 )/sqrt(2*length(dat(1).time));
             %error = sqrt( mean(dat(in(n).line).velU(:,doubleplot(1,:))).^2 + mean(dat(in(n).line).velU(:,doubleplot(2,:))).^2  )/sqrt(length(dat(1).time));
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


    set(ax(7),'ylim',[-2,20]);
end

end
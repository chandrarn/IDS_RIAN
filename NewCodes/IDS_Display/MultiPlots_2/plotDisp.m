%% Plot Displacement for Multiplots_2
function [saveDat] = plotDisp(dataDispl,ax,h,i,n,in,saveDat,dispSupress,...
    doubleplot,RMS,dat,data,pRel,CutPow,includeTemp,plotError,lnwdth)

% Loop through Arrays
for i = 1:1+in(n).doubleplot


%% We are doing the FFT    
% Remove data where the sine fit is invalid
for j=1:length(dataDispl(:,i));try dataDispl(j.*(pRel(j,i)<CutPow),i)=NaN;end;end 

% Plotting with errorbars?
if ~plotError
    saveDat(n).Displacement(:,i) = dataDispl(:,i);
    plot(ax(6),dat(1).impacts(1:size(data,2)),dataDispl(:,i),'-*','color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
else
    error = sqrt( nanmean(dat(in(n).line).velU(:,doubleplot(i,:)).^2).*(1e5/(2*pi*14500) ).^2 +...
        squeeze(RMS(i,:,n,100).*(1e5/(2*pi*14500))).^2)/sqrt(length(dat(1).time)) ; % convert error to centemeters displacement
    if any(dispSupress(n,:,i))
        dataDispl(dispSupress(n,:,i),i)=NaN;
    end
    saveDat(n).Displacement(:,i) = dataDispl(:,i);
    saveDat(n).DisplacementError(:,i) = error;
    saveDat(n).LMVelError(:,i) = ( nanmean(dat(in(n).line).velU(:,doubleplot(i,:))));
    errorbar(ax(6),dat(1).impacts(1:size(data,2)),dataDispl(:,i),error,'-*','color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
    saveDat(n).Displacement(1,i)
end

% Figure stuff
 ylabel(ax(6),'[cm]');set(ax(6),'ylim',[0,8]);
set(ax(6),'xlim',xlim);
if ~includeTemp
    xlabel(ax(6),'Impacts [cm]');
else
    xlabel(ax(17),'Impacts [cm]');
end

saveDat(n).Displacement(1,i);

end

end
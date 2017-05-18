%% Plot Phases For Multiplots_2
function [saveDat,phaseH] = plotPhase(dat,dataPhase,phaseSupress,plotError,...
    ax,h,n,data,CutPow,pRel,injParam,xlim,SigDev)

for i=1:size(doubleplot,1)
    if Analysis==1
        plot(dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,'-*','color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style);
        ylabel('Phase [deg]');
    elseif Analysis==2 && in(1).shot > 8129499
        for j=1:length(dataPhase(:,i));try dataPhase(j.*(pRel(j,i)<CutPow),i)=NaN;end;end
    %                         ax8Data = dataPhase(:,i);
    %                         if plotType==1
    %                             xindex=dat(1).impacts(1:size(data,2));
    %                         elseif plotType==2
    %                             xindex=dat(1).impacts(1*(i==1) +(size(data,2)+1)*(i==2):size(data,2)*(i==1)+(end)*(i==2));
    %                         end
        if any(phaseSupress(n,:,i))
            dataPhase(phaseSupress(n,:,i),i)=NaN;
        end
        if ~plotError
            saveDat(n).Phase(:,i) = dataPhase(:,i).*180./pi;
            phaseH(:,i,n)=plot(ax8,dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,'-*','color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
        else
            saveDat(n).Phase(:,i) = dataPhase(:,i).*180./pi;
            saveDat(n).PhaseError(:,i) = SigDev(n,:,i).*180./pi;
            phaseH(:,i,n)=errorbar(ax8,dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,SigDev(n,:,i).*180./pi,'-*','color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
        end                            
        ylabel(ax(8),'[deg]');set(ax(8),'ylim',[-400,400]);set(ax(8),'xticklabel',[]);
        set(ax(8),'xlim',xlim);
        if n==1
            saveDat(n).injPhase = injParam(:,3).*180./pi;
            saveDat(n).injPhase(1) = saveDat(n).injPhase(1)-20;
        plot(ax(8),[-60 60]',[injParam(1,3) injParam(1,3)]'.*180./pi-20,'--','LineWidth',lnwdth,'color',[0    0.4470    0.7410]);
        plot(ax(8),[-60 60]',[injParam(2,3) injParam(2,3)]'.*180./pi,'--','LineWidth',lnwdth,'color',[0.8500    0.3250    0.0980]);
        plot(ax(8),[-60 60]',[injParam(3,3) injParam(3,3)]'.*180./pi,'--','LineWidth',lnwdth,'color',[0.9290    0.6940    0.1250]);
        if in(1).shot >= 160728011 || in(1).shot == 151217024
            plot(ax(8),[-60 60]',[injParam(3,3) injParam(3,3)]'.*180./pi + 360,'--','LineWidth',lnwdth,'color',[0.9290    0.6940    0.1250]);
        end
        end

    elseif Analysis==2 && in(1).shot <= 8129499
        for j=1:length(dataPhase(:,1));try dataPhase(j.*(pRel(j,i)<CutPow),1)=NaN;end;end
        phaseH(1,1,n)=plot(ax8,dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
        ylabel(ax(8),'[deg]');set(ax8,'ylim',[0,400]);set(ax8,'xticklabel',[]);
        set(ax(8),'xlim',xlim);
        if n==1
        plot(ax(8),[0 60]',[injParam(1,3) injParam(1,3)]'.*180./pi,'--','LineWidth',lnwdth,'color',[    0    0.4470    0.7410]);
        plot(ax(8),[0 60]',[injParam(2,3) injParam(2,3)]'.*180./pi,'--','LineWidth',lnwdth,'color',[    0.8500    0.3250    0.0980]);
        end

    end
end

end
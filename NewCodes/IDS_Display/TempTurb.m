% Temp vs Turb

shot = 151217026;
tmRng = 50:200;
chnRng = 52:56;

load(['dat' num2str(shot) '10.mat']);

% find mean CII temp for each data point


h = figure('Name','CIII vs OII ScatterPlot'); hold on;

chSt = find(dat(1).peaks==chnRng(1),1);
chNd = find(dat(1).peaks==chnRng(end),1);

counter = 1;
for tm = tmRng
    for ch = chSt:chNd
        if ~isnan(dat(3).temp(tm,ch))
            TCIII(counter) = (dat(1).temp(tm,ch).*.7 + dat(3).temp(tm,ch).*.3);
        else
            TCIII(counter) = dat(1).temp(tm,ch);
        end
        TOII(counter) = dat(2).temp(tm,ch);
        if ~isnan(TCIII(counter)*TOII(counter))
            plot(TCIII,TOII,'*','MarkerSize',.5);
            counter = counter +1;
        end
    end
end


plot([12,48],[12,48],'-*');
plot([12,48],[12,48].*1.3,'k-*');
fun = inline('x(1)*(xdata)+0','x','xdata');
x =lsqcurvefit(fun,[1.2],TCIII,TOII)
plot([12,48],[12,48].*x,'g-*');
set(gca,'xlim',[10,50]);
set(gca,'ylim',[10,50]);
ylabel('OII');
xlabel('CIII');
title(['CIII vs OII scatter, Channels: ' num2str(chnRng(1)) '-' num2str(chnRng(end))]);

%{ 
alpha = mean(dat(1).param.PIX_SP).^2 * dat(1).param.c.^2 /( dat(1).param.kBoltz*11605)
dat(1).param.IonMass(1)/dat(1).param.LineLam(1).^2
dat(1).param.IonMass(2)/dat(1).param.LineLam(2).^2
slope ~ sigY(7.7)-8.1/(sigY(5.8)-6) asymptotically approaches 1.3, valid
for our sigY regime
%}


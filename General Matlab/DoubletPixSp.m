% calculate PIX_SP based on known doublet
%shot = 151217026;
shot = 160519019
timePoints = 75:250;
try dat;end
if ~strcmp(dat(1).title,'Shot 151217026');
load(['dat' num2str(shot) '10.mat']);
end

 Lambda = dat(1).param.LineLam(3)-dat(1).param.LineLam(1);
 Y1=dat(1).fit_par(timePoints,:,3);
Y3=dat(3).fit_par(timePoints,:,3);
Y1=averageNans(Y1);
Y3=averageNans(Y3);
Dist = Y3-Y1;
figure;
for i = 1:size(Y1,2)
    PIX_SP(i) = mean(Lambda./Dist(:,i));
    PIX_STD(i) = std(Lambda./Dist(:,i));
    errorbar(dat(1).peaks(i),PIX_SP(i),PIX_STD(i),'*');hold on;
end
plot(dat(1).peaks,PIX_SP)
set(gca,'ylim',[1.11,1.17].*1e-11);
xlabel('Channel Number');
ylabel('PIX_SP');
title('Doublet Callibrated PIX_SP');

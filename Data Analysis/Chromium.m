% We have chromium

import MDSplus.*
% load data 
shot = 150625140;
pixShot = 150401030; % if the chromium shot doesnt have pixSP
trimChan = 110:210;
trimTime = 5:12;
centerpoint = 22; % location of O II

%wavelengths
O = [ 464.913,465.084,464.181];
C = [464.742, 465.026, 465.147];
%Cr = [464.147];%,464.166,464.175,464.199,464.887,464.943];
Cu = [464.258,464.921];
N = [464.308];
W = [464.176,464.255];

addpath('T:\PhantomMovies');
data = importdata(['shot' int2str(shot) '.mat']);
% data = squeeze((mean(data(trimTime,:,trimChan),1)));
% data = squeeze(mean(data,2));
data = squeeze(sum(data(1,:,:),3));

hitTree = Tree('analysis3',pixShot);
%PIX = mean(NATIVEvalue(hitTree.getNode('\IDS_PIX_SP').getData()))*1e9;
PIX = 0.0114;
%PIX = 0.01166;
PIX = .0125
upperBound = 464.913+ PIX*(centerpoint-1);
lowerBound = 464.913 - PIX*(96-centerpoint);
x =(linspace(lowerBound,upperBound, 96));

for i = 1:96
    label{i} = num2str(x(i));
end


%plotting
figure;t(1,1)=plot(x,data(end:-1:1),'-*','LineWidth',3);
% set(gca,'xtick',[1:10:96]);
% set(gca,'xticklabel',label(1:10:96));
hold on;
set(gca,'FontSize',13);
title(['Spectral Lines for Shot: ' num2str(shot)]);
set(gca,'ytick',[]);
ylabel('Arb [ ]');
xlabel('Wavelength [nm]');

% make lines

% C = ( (C-O(1))./PIX + 75);
% Cr = ( (Cr-O(1))./PIX + 75);
% Cu = ( (Cu-O(1))./PIX + 75);
% N = ( (N-O(1))./PIX + 75);
% W = ( (W-O(1))./PIX + 96-21 );
% O = ( (O-O(1))./PIX + 96-21 );

for i = 1:length(O)
    t(2,i)=plot([O(i),O(i)],[.8e5,1.1e5],'c','LineWidth',3);
end
for i = 1:length(C)
    t(3,i)=plot([C(i),C(i)],[.8e5,1.1e5],'k','LineWidth',3);
end
% for i = 1:length(Cr)
%     t(4,i)=plot([Cr(i),Cr(i)],[65,90],'y','LineWidth',3);
% end
for i = 1:length(Cu)
    t(5,i)=plot([Cu(i),Cu(i)],[.8e5,1.1e5],'r','LineWidth',3);
end
for i = 1:length(N)
    t(6,i)=plot([N(i),N(i)],[.8e5,1.1e5],'g','LineWidth',3);
end
for i = 1:length(W)
    t(7,i)=plot([W(i),W(i)],[.8e5,1.1e5],'m','LineWidth',3);
end

legend([t(1,1),t(2,1),t(3,1),t(5,1),t(6,1),t(7,1)],{'Spectra','O II','C III','Cr I','Cu I','N II','W I'})

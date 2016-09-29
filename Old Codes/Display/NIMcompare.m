% plotting routine for NIMROD/MC_IDS comparison

Cdats = 150303021; % CIDS shot
Mdats = 150310023; % MCIDS shot
Csft = 0;
Msft = 2;
Crng = 4:6; % CIDS framerange 
Mrng = 1:2;% MCIDS framerange
d = [12.12047564931123D-12,12.1220456743384D-12,12.12658131000742D-12,12.12390593011439D-12,12.11947900981399D-12,12.1209830914616D-12,12.12409138273652D-12,12.10731780684983D-12,12.11167308183358D-12,12.11762225997158D-12,12.11986463380678D-12,12.11751362646397D-12,12.1216173840457D-12,12.12270495324322D-12,12.1162896975129D-12,12.11346808951918D-12,12.1304716531523D-12,12.12288660358735D-12,12.11166676472455D-12,12.10330852924166D-12,12.10251235118542D-12,12.13639099603379D-12,12.09421711208199D-12,12.10329555685615D-12,12.09323541378095D-12,12.0905664608748D-12,12.09004697806257D-12];
% new PIX_SP value, use for both C and M
addpath('T:\IDS\Data Repository');
addpath('E:\Cines');
Cdat = importdata(['dat' num2str(Cdats) '10.mat']);
Mdat = importdata(['dat' num2str(Mdats) '10.mat']);

Mfac = d./Mdat.param.PIX_SP'; % MCIDS correction factor
[P,S] = polyfit(Mdat.impacts,d',1);
d1 = polyval(P,Cdat.impacts-13,S); 
Cfac = d1./Cdat.param.PIX_SP; % CIDS correction factor
% 
% for i = Crng
%     Cdat.vel(i,:) = Cdat.vel(i,:) .* Cfac';
% end
% for i = Mrng
%     Mdat.vel(i,:) = Mdat.vel(i,:).*Mfac;
% end

% SURF THE VEL
figure;

subplot(2,1,1); hold on;
surf(Mdat.vel+Msft); shading interp; view([ 0 90]);
title(['MCIDS Shot: ' num2str(Mdats)],'fontsize',13);
set(gca,'fontsize',13);
set(gca,'xlim',[1,27]);
set(gca,'ylim',[1,length(Mdat.time)]);
xlabel('Channel','fontsize',13);
ylabel('Frame','fontsize',13);
subplot(2,1,2); hold on;
title(['CIDS Shot: ' num2str(Cdats)],'fontsize',13);
surf(Cdat.vel+Csft); shading interp; view([ 0 90]);
set(gca,'fontsize',13);
set(gca,'xlim',[1,168]);
set(gca,'ylim',[1,length(Cdat.time)]);
xlabel('Channel','fontsize',13);
ylabel('Frame','fontsize',13);

% matlab color vector
rgb= [   0         0    1.0000
         0    0.5000         0
    1.0000         0         0
         0    0.7500    0.7500
    0.7500         0    0.7500
    0.7500    0.7500         0
    0.2500    0.2500    0.2500];
figure;

% velocity
h(1)=subplot(4,1,1); hold on;
plot(Cdat.impacts-13,Cdat.vel(Crng,:)+Csft);
plot(Mdat.impacts, Mdat.vel(Mrng,:)+Msft,'-*');
ylabel('Velocity, [Km/S]','FontSize',13);


% temperature
h(2)=subplot(4,1,2); hold on;
plot(Cdat.impacts-13,Cdat.temp(Crng,:));
plot(Mdat.impacts, Mdat.temp(Mrng,:),'-*');
ylabel('Temperature, [Ev]','FontSize',13);

% Rel Int
h(3)=subplot(4,1,3); hold on;
plot(Cdat.impacts-13,Cdat.int(Crng,:));
plot(Mdat.impacts, Mdat.int(Mrng,:),'-*');
ylabel('Intensity, [arb]','FontSize',13);
xlabel('Impact Parameters','FontSize',13);

% Current
h(4)=subplot(4,1,4); hold on;
plot(Cdat.iinjaTime,Cdat.iinja,Cdat.ItorTime,Cdat.Itor,Mdat.ItorTime,Mdat.Itor);
for i = Crng
    plot([Cdat.time(i)*1e-6;Cdat.time(i)*1e-6],[-10,10],'color',rgb(i-Crng(1)+1,:),'LineWidth',3);
end
for i = Mrng
    plot([Mdat.time(i)*1e-6;Mdat.time(i)*1e-6],[-10,10],'--*','color',rgb(i-Mrng(1)+1,:),'LineWidth',3);
end
xlim([.5,2.5]);
ylim([-20, 30]);
xlabel('Time, [S]','FontSize',13);
ylabel('Current, [kV]','FontSize',13);

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, .985,['\bf MCIDS Shot: ' num2str(Mdats) ' "*", CIDS Shot: ' num2str(Cdats) ' "-"'],'FontSize', 13, 'HorizontalAlignment','center','VerticalAlignment', 'top');
for i = 1:4
set(h(i),'FontSize',13);
end
set(h(1),'xtick',[]);
set(h(2),'xtick',[]);
set(h(1),'Position',[.1,.75,.85,.203]);
set(h(2), 'Position',[.1,.523,.85,.203]);
set(h(3),'Position',[.1,.29,.85,.203]);
set(h(4),'Position',[.1,.09,.85,.1359]);





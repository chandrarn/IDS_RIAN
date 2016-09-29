% plotting routine for C_IDS/MC_IDS comparison

close all;
dats1 = 150401016; % CIDS shot
dats2 = 150401018; % MCIDS shot
sft1 = 0;%2
sft2 = 0;
rng1 = 6:11; % CIDS framerange 
rng2 = 6:11;% MCIDS framerange
chanrng1 = [7:27];[1:23];
chanrng2 = [7:27];[1:23];
chanrT1 = [1:length(chanrng1)];-chanrng1(1)+1;;%[6:26]-chanrng1(1)+1;%[ 1+30:168-chanrng1(1)];
chanrT2 = [4:length(chanrng2)];-chanrng2(1)+1;;%[9:27]-chanrng2(1)+1;%[ 1+4:28-chanrng2(1)]; % special range for temp, frequently noisy

%d = [12.12047564931123D-12,12.1220456743384D-12,12.12658131000742D-12,12.12390593011439D-12,12.11947900981399D-12,12.1209830914616D-12,12.12409138273652D-12,12.10731780684983D-12,12.11167308183358D-12,12.11762225997158D-12,12.11986463380678D-12,12.11751362646397D-12,12.1216173840457D-12,12.12270495324322D-12,12.1162896975129D-12,12.11346808951918D-12,12.1304716531523D-12,12.12288660358735D-12,12.11166676472455D-12,12.10330852924166D-12,12.10251235118542D-12,12.13639099603379D-12,12.09421711208199D-12,12.10329555685615D-12,12.09323541378095D-12,12.0905664608748D-12,12.09004697806257D-12];
% new PIX_SP value, use for both C and M
%surpress Temperature, frequently noisy
surpT1 = 0;
surpT2 = 0;
addpath('T:\IDS\Data Repository');
addpath('E:\Cines');
dat1 = importdata(['dat' num2str(dats1) '10.mat']);
dat2 = importdata(['dat' num2str(dats2) '10.mat']);

%flip CIDS impacts
%dat1.impacts = dat1.impacts(end:-1:1);
% trim range if important
if ~isempty(chanrng1)
    dat1=trimRange(dat1,chanrng1);
    dat2=trimRange(dat2,chanrng2);
end


if isempty(chanrT1)
    chanrT1 = 1:length(dat1.impacts);
end
if isempty(chanrT2)
    chanrT2 = 1:length(dat2.impacts);
end
% for when you're REALLY done:
 %dat2.temp(11,8) = NaN;
 %dat2.vel(7,1:3) = NaN;
 %dat1.temp(1,3:5)= nan;
 %dat2.vel(7,2) = nan;
 %dat1.time(3) = 2*dat1.time(2) -dat1.time(1);
% Mfac = d./dat1.param.PIX_SP'; % MCIDS correction factor
% [P,S] = polyfit(dat1.impacts,d',1);
% d1 = polyval(P,dat2.impacts-13,S); 
% Cfac = d1./dat2.param.PIX_SP; % CIDS correction factor
% 
% for i = rng1
%     dat2.vel(i,:) = dat2.vel(i,:) .* Cfac';
% end
% for i = rng2
%     dat1.vel(i,:) = dat1.vel(i,:).*Mfac;
% end
%Check if IDS or CIDS
if dats1 <150310021
    dat1T = 'CIDS';
    impshft1 = -13;
else
    dat1T = 'MCIDS';
    impshft1 = 0;
end
if dats2 <150310021
    dat2T = 'CIDS';
    impshft2 = -13;
else
    dat2T = 'MCIDS';
    impshft2 = 0;
end

% SURF THE VEL
figure;

subplot(2,1,1); hold on;
surf(dat1.vel+sft1); shading interp; view([ 0 90]); box on;

title([dat1T ' Shot: ' num2str(dats1)],'fontsize',13);
set(gca,'fontsize',13);
set(gca,'xlim',[1,length(dat1.peaks)]);
set(gca,'ylim',[1,length(dat1.time)]);
xlabel('Channel','fontsize',13);
ylabel('Frame','fontsize',13);
subplot(2,1,2); hold on;
title([dat2T ' Shot: ' num2str(dats2)],'fontsize',13);
surf(dat2.vel+sft2); shading interp; view([ 0 90]);
set(gca,'fontsize',13);
set(gca,'xlim',[1,length(dat2.peaks)]);
set(gca,'ylim',[1,length(dat2.time)]);
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
plot(dat1.impacts+impshft1,dat1.vel(rng1,:)+sft1);
plot(dat2.impacts+impshft2, dat2.vel(rng2,:)+sft2,'-*');
ylabel('Velocity [km/s]','FontSize',13);
ylim([-20 20]);
%title([num2str(dats1) ': "-", ' num2str(dats2) ': "-*-"'],'fontsize',13);
set(h(1),'Xlim',[-30 50]); 
grid on;
box on;

% temperature
h(2)=subplot(4,1,2); hold on;
if surpT1 == 0
    plot(dat1.impacts(chanrT1)+impshft1,dat1.temp(rng1,chanrT1));
else
    text(0,20,[dat1T ' SURPRESSED'],'fontweight','bold');
end
if surpT2 == 0
    plot(dat2.impacts(chanrT2)+impshft2, dat2.temp(rng2,chanrT2),'-*');
else
    text(20,60,[dat2T ' SURPRESSED']);
end
% if they're both surpressed
if surpT1 == 1 && surpT2 == 1
    %plot(dat2.impacts+impshft2, dat2.temp(rng2,:),'-*');
    set(h(2),'ytick',[]);
    text(.23,.5,'\bf TEMPERATURE SURPRESSED','fontsize',13);
end
ylabel('Temperature [eV]','FontSize',13);
set(h(2),'Xlim',[-30 50]); % maintain common plotting limits
grid on;
box on;

% Rel Int
h(3)=subplot(4,1,3); hold on;
plot(dat1.impacts+impshft1,dat1.int(rng1,:));
plot(dat2.impacts+impshft2, dat2.int(rng2,:),'-*');
ylabel('Intensity [arb]','FontSize',13);
xlabel('Impact Parameters','FontSize',13);
set(h(3),'Xlim',[-30 50]); 
grid on
box on;

% Current
h(4)=subplot(4,1,4); hold on;

for i = rng1
    plot([dat1.time(i)*1e-6;dat1.time(i)*1e-6],[-10,10],'color',rgb(i-rng1(1)+1,:),'LineWidth',3);
end
for i = rng2
    plot([dat2.time(i)*1e-6;dat2.time(i)*1e-6],[-10,10],'--*','color',rgb(i-rng2(1)+1,:),'LineWidth',3);
end
plot(dat2.iinjaTime,dat2.iinja)
d2 = plot(dat2.ItorTime,dat2.Itor,'color',[0 .5 0]);
d1 = plot(dat1.ItorTime,dat1.Itor,'r');
l=legend([d1,d2],'-','-*-');
set(h(4),'box','on')
set(l,'fontsize',9)
set(l,'fontweight','bold')
xlim([.5,2.5]);
ylim([-30, 30]);
xlabel('Time, [ms]','FontSize',13);
ylabel('Current [kA]','FontSize',13);





ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, .985,['\bf ' dat1T ' Shot: ' num2str(dats1) ' "-", ' dat2T ' Shot: ' num2str(dats2) ' "-*-"'],'FontSize', 13, 'HorizontalAlignment','center','VerticalAlignment', 'top');
for i = 1:4 
set(h(i),'FontSize',13);
end
set(h(1),'xticklabel',[]);
set(h(2),'xticklabel',[]);
set(h(1),'Position',[.1,.75,.85,.203]);
set(h(2), 'Position',[.1,.523,.85,.203]);
set(h(3),'Position',[.1,.29,.85,.203]);
set(h(4),'Position',[.1,.09,.85,.1359]);
set(h(1),'linewidth',2)
set(h(2),'linewidth',2)
set(h(3),'linewidth',2)





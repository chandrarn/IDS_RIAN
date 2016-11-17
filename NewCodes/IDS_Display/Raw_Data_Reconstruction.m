% Plot figure to show fit, error, etc. 
% using 0-120-240, impact 32.8, index 9

shot= 160728011;
 addpath('T:\IDS\Data Repository');
 input=load(['dat' num2str(shot) '10.mat']); % Real HIT-SI Data
 c=colormap;
 n=1;
i=9;
velShift=0;
 timebound = [1.65,2.00];
    chan_ranget = [9:24];
    chan_rangep = [45:60];
    chan_range = [chan_ranget, chan_rangep];
    timeScale = 1e-3;
dat = trimRange(input.dat, chan_range, 0,timebound.*(1./timeScale),[]);

 % Find where each fiber bundle begins and ends.
  dat(in(n).line).vel = averageNans(dat(in(n).line).vel)+velShift; % remove nans
  
  data(1:length(dat(1).time),:) = dat(in(n).line).vel(:,1:(length(dat(1).impacts))/2)+velShift;
                data(length(dat(1).time)+1:2*length(dat(1).time),:) = ...
                   dat(in(n).line).vel(:,(length(dat(1).impacts)/2)+1:end)+velShift;
      
doubleplot(1,:) = 1:(length(dat(1).impacts))/2;
doubleplot(2,:) = (length(dat(1).impacts)/2)+1:length(dat(1).impacts);

% CALC UPPER FIBER

pRel = zeros(length(doubleplot),2);
dPar = zeros(length(in),size(doubleplot,2),size(doubleplot,2),4);

signal1 = dat(in(n).line).vel(:,doubleplot(1,i));
Fsamp = 1/(mean(diff(dat(1).time.*(timeScale.*1e-3))));
offset = nanmean(signal1);
amp = max(signal1)-offset;
freq = 14500;
phase = pi/2;
guess(i,:,n,1) = [offset,amp,phase,freq];
[param(i,2:5,n),data(1:length(dat(1).time),i)] = sine_fit( ...
    dat(1).time'.*(timeScale.*1e-3),signal1',[nan,nan,nan,freq], ...
    [offset,amp,phase,freq],0);
if param(i,3,n)<0 % 180degree phase
    param(i,3,n)=-param(i,3,n);
    param(i,4,n)=param(i,4,n)+pi;
    disp([' WARNING: NEGATIVE AMPLITUDE @ n=' num2str(n) ', Line 1, Impact: ' num2str(dat(1).impacts(i))]);
end

% CALC LOWER FIBER

pRel = zeros(length(doubleplot),2);
dPar = zeros(length(in),size(doubleplot,2),size(doubleplot,2),4);

signal = dat(in(n).line).vel(:,doubleplot(2,i));
Fsamp = 1/(mean(diff(dat(1).time.*(timeScale.*1e-3))));
offset = nanmean(signal);
amp = max(signal)-offset;
freq = 14500;
phase = pi/2;
guess(i,:,n,2) = [offset,amp,phase,freq];
[param(i,7:10,n),data(1:length(dat(1).time),i)] = sine_fit( ...
    dat(1).time'.*(timeScale.*1e-3),signal',[nan,nan,nan,freq], ...
    [offset,amp,phase,freq],0);
if param(i,8,n)<0 % 180degree phase
    param(i,8,n)=-param(i,8,n);
    param(i,9,n)=param(i,9,n)+pi;
    disp([' WARNING: NEGATIVE AMPLITUDE @ n=' num2str(n) ', Line 1, Impact: ' num2str(dat(1).impacts(i))]);
end


figure; plot(dat(1).time'.*(timeScale.*1e-3),signal1,dat(1).time'.*(timeScale.*1e-3),signal,'linewidth',2);
hold on;
time = linspace(dat(1).time(1).*(timeScale.*1e-3),dat(1).time(end).*(timeScale.*1e-3),length(dat(1).time).*10);
plot(time,SineFitLM(time,param(i,2:5,n)),'--','color',colorOrder(1,:),'linewidth',2);
plot(time,SineFitLM(time,param(i,7:10,n)),'--','color',colorOrder(2,:),'linewidth',2);
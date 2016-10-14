%% Convert NIMROD postprocessing output to normal IDS form

Nimshot = 160609202;
line = [1];
hitsi3=1;
useTree=1;

%% Param Stuff
if ~exist('nimsave','var')
    load(['T:\IDS\Data Repository\NIMROD\nim' num2str(Nimshot) '.mat']);
end
shot = str2num(nimsave.nimIn.wform_fname(2:10));
% temp override
shot = 160622021;
[dat(1).param, options] = loadParams(shot, line, hitsi3, useTree);
shot = str2num(nimsave.nimIn.wform_fname(2:10));
impacts=load('T:\RChandra\NewCodes\Geometry\impacts5.mat');
dat(1).param.impacts=impacts.impacts;

dat(1).param=param;
%% Injector Stuff
dat(1).iinja = 1e-3 * nimsave.inj.Ia';
dat(1).iinjb = 1e-3 * nimsave.inj.Ib';
dat(1).iinjc = 1e-3 * nimsave.inj.Ic';
dat(1).iinjaTime = 1e3*nimsave.time';
dat(1).iinjbTime = 1e3*nimsave.time';
dat(1).iinjcTime = 1e3*nimsave.time';
Itor=nimsave.Bprob.surf.Itor_spa_avg;
ItorTime = nimsave.time;
dat(1).ItorTime = 1e3 * ItorTime;
dat(1).Itor = 1e-3 * Itor;


%% Data
dat(1).vel = 1e-3*nimsave.IDS.vIDS(1:72,:)';
dat(1).temp = nimsave.IDS.TIDS(1:72,:)';
dat(1).time=1e6*nimsave.time'; %Time in nS
dat(1).int=zeros(size(dat(1).vel));
dat(1).fit_par=zeros([size(dat(1).vel) 6]);
dat(1).bounds=zeros([size(dat(1).vel) 4]);
dat(1).guesses=zeros([size(dat(1).vel) 6]);
dat(1).int=zeros([size(dat(1).vel) 33]);
dat(1).peaks=1:72;


dat(1).impacts = impacts.impacts;
dat(1).shotRef=shot;
dat(1).label1= 'Impact Parameter [cm]';
dat(1).label2='Major Radius [cm]';
dat(1).title=['Shot 8' num2str(shot)];



%% Saving
dat=dat(1);
saveName = ['T:\IDS\Data Repository\dat8' num2str(shot) '10.mat']
save(['' saveName],'dat');

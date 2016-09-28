%% this program compares the IMP data to Taylor
%% Updated to new MDSplus by Rian Chandra

clear all; close all;
import MDSplus.*
addpath('U:\bvictor\Matlab\filters');% accesses over the network

% save plots?
%newdir = ['S:\Matlab\high_freq_pop\figures\parameters'];
%mkdir(newdir); % DONT OVERWITE BRIAN'S STUFF ON ACCIDENT


%% global variables and inputs are here
% these constants scale the Taylor data
mu0 = 4*pi*1e-7;
kboltz = 1.38e-23;

f_order = 4; % order of all the filters
low_filt = 23000; % low pass frequency
high_filt = 27000; % high pass

% put this option in for probes, too
filterdata = 0; % 2 = bandpass, 1 = lowpass, 0 = no filter
% 3 = lowpass at half the Nyquist
% 4 = removes a frequency range from the data

% r is in the i direction
% p is in the j direction
% t is in the k direction
dir = ['R'; 'P'; 'T'];

shots = [129495, 129496, 129793,   129440, 129441, 129443, 129446, 129449, 129450, 129451]; % input('What shot? ');
%[a, status] = mdsopen('landau.hit::hitsi', shot);
HitTree = Connection('landau.hit'); % New MDS version for thin client connection
%HitTree.openTree('hitsi',shot);


% define variables
ts = -0.0001; % start time for new time base
ts2 = ts - 0.0002; % start time for cropping data in
te = 0.005; % end time for new time base
te2 = te + 0.0002; % end time for cropping data in
dtn = 0.4e-6; % dt for new time base (2.5 MHz)
tbase = ts: dtn: te; % new time base
color = ['-r'; '-g'; '-b'; '-c'; '-m'; '-k';
    ':r'; ':g'; ':b'; ':c'; ':m'; ':k'];

%----------------------------------------------------------------%
%% this brings in the data IN A LOOP
for i = 1:length(shots)
    HitTree.openTree('hitsi',shots(i));% open shot
    
    [shot(i).t.vxinj, shot(i).dt.vxinj, shot(i).sig.vxinj] = gen_data_in('\v_inj_x',HitTree);
    [shot(i).t.vyinj, shot(i).dt.vyinj, shot(i).sig.vyinj] = gen_data_in('\v_inj_y',HitTree);
    [shot(i).t.psixinj, shot(i).dt.psixinj, shot(i).sig.psixinj] = gen_data_in('\psi_inj_x',HitTree);
    [shot(i).t.psiyinj, shot(i).dt.psiyinj, shot(i).sig.psiyinj] = gen_data_in('\psi_inj_y',HitTree);
    [shot(i).t.ixinj, shot(i).dt.ixinj, shot(i).sig.ixinj] = gen_data_in('\i_inj_x',HitTree);
    [shot(i).t.iyinj, shot(i).dt.iyinj, shot(i).sig.iyinj] = gen_data_in('\i_inj_y',HitTree);
    [shot(i).t.itor, shot(i).dt.itor, shot(i).sig.itor] = gen_data_in('\i_tor_spaavg',HitTree);
    [shot(i).t.den, shot(i).dt.den, shot(i).sig.den] = gen_data_in_smart('\n_avg_s1',HitTree);

    [shot(i).t.itorsm, shot(i).dt.itorsm, shot(i).sig.itorsm] = gen_data_in('sihi_smooth(\i_tor_spaavg)',HitTree);
    [shot(i).t_qinj, shot(i).dt_qinj, shot(i).sig_qinj] = ...
        gen_data_in('sihi_smooth(sigadd_quad(\i_inj_x,\i_inj_y))',HitTree);
    shot(i).int_itor = interp1(shot(i).t.itorsm, shot(i).sig.itorsm, tbase);
    shot(i).int_iquad = interp1(shot(i).t_qinj, shot(i).sig_qinj, tbase);
    shot(i).i_ratio = abs(shot(i).int_itor./shot(i).int_iquad);

    % calculate the signal amplitudes here
    shot(i).hil_ixinj = hilbert_bsv(shot(i).sig.ixinj');
    shot(i).hil_iyinj = hilbert_bsv(shot(i).sig.iyinj');
    shot(i).hil_psixinj = hilbert_bsv(shot(i).sig.psixinj');
    shot(i).hil_psiyinj = hilbert_bsv(shot(i).sig.psiyinj');
    shot(i).hil_vxinj = hilbert_bsv(shot(i).sig.vxinj');
    shot(i).hil_vyinj = hilbert_bsv(shot(i).sig.vyinj');

    shot(i).mag_ixinj = sqrt(shot(i).sig.ixinj.^2 + shot(i).hil_ixinj'.^2);
    shot(i).mag_iyinj = sqrt(shot(i).sig.iyinj.^2 + shot(i).hil_iyinj'.^2);
    shot(i).mag_psixinj = sqrt(shot(i).sig.psixinj.^2 + shot(i).hil_psixinj'.^2);
    shot(i).mag_psiyinj = sqrt(shot(i).sig.psiyinj.^2 + shot(i).hil_psiyinj'.^2);
    shot(i).mag_vxinj = sqrt(shot(i).sig.vxinj.^2 + shot(i).hil_vxinj'.^2);
    shot(i).mag_vyinj = sqrt(shot(i).sig.vyinj.^2 + shot(i).hil_vyinj'.^2);

    % smooth the signals over an injector cycle
    % Sometimes, mdsvalue works. Sometimes, it doesnt. 
    % If this cycles too many times, close matlab, put on coat, go to
    % shultzy's until caring about the problem is eliminated. 
    mdsvalue('_ixinj = make_signal($1, *, $2)', shot(i).mag_ixinj, shot(i).t.ixinj);
    mdsvalue('_iyinj = make_signal($1, *, $2)', shot(i).mag_iyinj, shot(i).t.iyinj);
    mdsvalue('_psixinj = make_signal($1, *, $2)', shot(i).mag_psixinj, shot(i).t.psixinj);
    mdsvalue('_psiyinj = make_signal($1, *, $2)', shot(i).mag_psiyinj, shot(i).t.psiyinj);
    mdsvalue('_vxinj = make_signal($1, *, $2)', shot(i).mag_vxinj, shot(i).t.vxinj);
    mdsvalue('_vyinj = make_signal($1, *, $2)', shot(i).mag_vyinj, shot(i).t.vyinj);

    % sihi_ixinj = mdsvalue('data(sihi_smooth(_ixinj))');
    % sihi_iyinj = mdsvalue('data(sihi_smooth(_iyinj))');
    % sihi_psixinj = mdsvalue('data(sihi_smooth(_psixinj))');
    % sihi_psiyinj = mdsvalue('data(sihi_smooth(_psiyinj))');
    % sihi_vxinj = mdsvalue('data(sihi_smooth(_vxinj))');
    % sihi_vyinj = mdsvalue('data(sihi_smooth(_vyinj))');

    shot(i).aihi_ixinj = shot(i).mag_ixinj;
    shot(i).sihi_iyinj = shot(i).mag_iyinj;
    shot(i).sihi_psixinj = shot(i).mag_psixinj;
    shot(i).sihi_psiyinj = shot(i).mag_psiyinj;
    shot(i).sihi_vxinj = shot(i).mag_vxinj;
    shot(i).sihi_vyinj = shot(i).mag_vyinj;
end
mdsclose;
%% plotting routines here

ts_den = -.2;
% ts_den = .00092; %start plot from here b/c of fringe jumps
[aa, I_den] = min((shot(i).t.den - ts_den).^2);

% if shot == 129175
%     tps = -0.01;
%     tpe2 = 2.3;
% 
%     ytick_a2 = [0 .25 .5];
%     ytick_l2 = [-.1 .6];
%     
% %     ytick_a3 = [0 4 8];
% %     ytick_l3 = [-1 10];
%     ytick_a3 = [0 15 30];
%     ytick_l3 = [-5 31];
%     
%     ytick_a4 = [0 .9 1.8];
%     ytick_l4 = [-.2 2];
%     ytick_a5 = [0 1.5 3];
%     ytick_l5 = [-1 4];
%     s_ratio = .000683;
%     e_ratio = .00183;
% 
% elseif shot == 128436
%     tps = -0.01;
%     tpe2 = 2.5;
% 
%     ytick_a2 = [0 .25 .5];
%     ytick_l2 = [-.1 .6];
%     
% %     ytick_a3 = [0 5 10];
% %     ytick_l3 = [-1 12];
%     ytick_a3 = [0 15 30];
%     ytick_l3 = [-5 31];
%     
%     ytick_a4 = [0 1.5 3];
%     ytick_l4 = [-.5 4];
%     ytick_a5 = [0 1.5 3];
%     ytick_l5 = [-1 4];
%     s_ratio = .00067;
%     e_ratio = .00198;
%     
% elseif shot == 126482
%     tps = -0.01;
%     tpe2 = 2.5;
% 
% %     ytick_a2 = [0 .5 1];
% %     ytick_l2 = [-.1 1.2];
%     ytick_a2 = [0 .5 1];
%     ytick_l2 = [-.1 1.2];
%     
%     ytick_a3 = [0 20 40];
%     ytick_l3 = [-8 43];
%     
%     ytick_a4 = [0 2 4];
%     ytick_l4 = [-1 5];
%     ytick_a5 = [0 1.5 3];
%     ytick_l5 = [-.5 4];
%     s_ratio = .000755;
%     e_ratio = .00198;
%     
% elseif shot == 122385
%     tps = -0.01;
%     tpe2 = 2.4;
% 
%     ytick_a2 = [0 .5 1];
%     ytick_l2 = [-.1 1.4];
%     
%     ytick_a3 = [0 25 50];
%     ytick_l3 = [-5 57];
% %     ytick_a3 = [0 10 20];
% %     ytick_l3 = [-1 24];
% 
%     ytick_a4 = [0 3 6];
%     ytick_l4 = [-1 8];
%     ytick_a5 = [0 1.5 3];
%     ytick_l5 = [-.5 4];
%     s_ratio = .000755;
%     e_ratio = .001753;
%     
% else % copy of above
    tps = -0.01;
    tpe2 = 2.4;
    
    ytick_a2 = [0 .5 1];
    ytick_l2 = [-.1 1.4];
    
    ytick_a3 = [0 25 50];
    ytick_l3 = [-5 57];
%     ytick_a3 = [0 10 20];
%     ytick_l3 = [-1 24];

    ytick_a4 = [0 3 6];
    ytick_l4 = [-1 8];
    ytick_a5 = [0 1.5 3];
    ytick_l5 = [-.5 4];
    s_ratio = .000755;
    e_ratio = .001826+.0005;
    
    
    currentLim=[-100 10];
    quadLim = [-5 40];
    ampLim= [0 4];
% end


x_ticks = [0.5 1 1.5 2];

scrsz = get(0,'ScreenSize');
% declair all the plots outside of the loop, plot the values inside it
% voltage flux and current
h1 = figure(8);     %[Left edge, Bottom edge, Width, Height]
set(h1, 'Position', [1 scrsz(4)/20 scrsz(3)/3 scrsz(4)/1.5],'Color',[1 1 1]);
set(h1, 'Name', 'params');
% set(h1,'itle','CURRENTS');
% title('Current Measurements');
hold('all');
clf;

axes1 = axes('Parent',h1,'YTick',[-80 -40 0],'XTickLabel','',...
    'XTick',x_ticks,...
    'Position',[0.15 0.75 0.84 0.17],...
    'FontSize',16);
title('Example Measurements (Negative Shots)');

box('on'); grid('on'); hold('all');
xlim([tps tpe2]);
ylim(currentLim);
ylabel('I_{tor} [kA]', 'FontSize', 16);
%legend('X-inj', 'Y-inj', 'Location', 'SouthWest', 'Orientation', 'Horizontal');
%text(1.51, -100, int2str(shots(1)), 'Color', 'k', 'Fontsize', 14);

axes2 = axes('Parent',h1,'YTick', [0 15 30],'XTickLabel','',...
    'XTick',x_ticks,...
    'Position',[0.15 0.575 0.84 0.17],...
    'FontSize',16);
box('on'); grid('on'); hold('all');
xlim([tps tpe2]);
ylim(quadLim);
ylabel('I_{inj} [kA]', 'FontSize', 16);


axes3 = axes('Parent',h1,'YTick',[0 1.5 3],...
    'Position',[0.15 0.40 0.84 0.17],...
    'FontSize',16);
box('on');
grid('on');
hold('all');
xlim([tps tpe2]);
ylim([0 3.5]);
ylabel('I_{tor}/I_{inj}', 'FontSize', 16);
xlabel('time [ms]', 'FontSize', 18);

% axes3 = axes('Parent',h1,'YTick',[0 1.5 3],'XTickLabel','',...
%     'XTick',x_ticks,...
%     'Position',[0.15 0.45 0.84 0.17],...
%     'FontSize',16);
% box('on'); grid('on'); hold('all');
% xlim([tps tpe2]);
% ylim(ampLim);
% ylabel('I_{tor}/I_{inj}', 'FontSize', 16);
% xlabel('time [ms]', 'FontSize', 18);

% Plot the plots
% temp1=shot(4).i_ratio
% size(shot(4).i_ratio)
% shot(3).i_ratio=-1.*(shot(3).i_ratio);
% shot(4).i_ratio=-1.*(shot(4).i_ratio);
% shot(5).i_ratio=(shot(5).i_ratio);
[aa, Is] = min((tbase - s_ratio).^2)
[aa, Ie] = min((tbase - e_ratio).^2)

for i =1:length(shots)
    i
    plot(shot(i).t.itor*1e3, shot(i).sig.itor*1e-3, 'LineWidth', 2,'Parent',axes1);

    plot(tbase*1e3, shot(i).int_iquad*1e-3, 'LineWidth', 2,'Parent',axes2);
    
    
%     TEMP=size(shot(i).i_ratio)
      plot(tbase(Is: Ie)*1e3, shot(i).i_ratio(Is: Ie), 'LineWidth', 2,'Parent',axes3);
%     TEMP=size(shot(i).i_ratio)
end



% axes3 = axes('Parent',h1,'YTick',[0 1.5 3],'XTickLabel','',...
%     'XTick',x_ticks,...
%     'Position',[0.15 0.45 0.84 0.17],...
%     'FontSize',16);
% box('on'); grid('on'); hold('all');
% size(shot)
% size(shot(1))
% size(tbase)
% size(tbase(Is: Ie))
% size(shot(1).i_ratio)
% size(shot(1).i_ratio(1,Is: Ie))
% plot(tbase(Is: Ie)*1e3, shot(1).i_ratio(1,Is: Ie), 'LineWidth', 2,'Parent',axes3);
% plot(tbase(Is: Ie)*1e3, shot(2).i_ratio(1,Is: Ie), 'LineWidth', 2,'Parent',axes3);
% plot(tbase(Is: Ie)*1e3, shot(3).i_ratio(1,Is: Ie), 'LineWidth', 2,'Parent',axes3);
% xlim([tps tpe2]);
% ylim(ampLim);
% ylabel('I_{tor}/I_{inj}', 'FontSize', 16);
% xlabel('time [ms]', 'FontSize', 18);


%  plot(tbase(Is: Ie)*1e3, shot(1).i_ratio(1,Is: Ie), 'LineWidth', 2,'Parent',axes3);
%  plot(tbase(Is: Ie)*1e3, shot(2).i_ratio(1,Is: Ie), 'LineWidth', 2,'Parent',axes3);
%  plot(tbase(Is: Ie)*1e3, shot(3).i_ratio(1,Is: Ie), 'LineWidth', 2,'Parent',axes3);
%plot(shot(1).t.itor*1e3, shot(1).sig.itor*1e-3, 'k', 'LineWidth', 2,'Parent',axes5);





h=legend(axes1,cellstr(num2str(shots')),'Location','SouthWest','Orientation','vertical');
set(h,'FontSize',10);
fig_save = getframe(h1);
[Xfig, mapfig] = frame2im(fig_save);
%imwrite(Xfig, [newdir '\params' int2str(shot) '.png']);
% 
% 
% % voltage flux and current
% h1 = figure(9);     %[Left edge, Bottom edge, Width, Height]
% set(h1, 'Position', [1 scrsz(4)/20 scrsz(3)/3 scrsz(4)/1.5],'Color',[1 1 1]);
% set(h1, 'Name', 'params');
% hold('all');
% clf;
% axes1 = axes('Parent',h1,'YTick',[0 250 500],'XTickLabel','',...
%     'XTick',x_ticks,...
%     'Position',[0.15 0.8 0.84 0.17],...
%     'FontSize',16);
% box('on');
% grid('on');
% hold('all');
% plot(t.vxinj*1e3, sihi_vxinj, 'r', 'LineWidth', 2,'Parent',axes1);
% plot(t.vyinj*1e3, sihi_vyinj, '--b', 'LineWidth', 2,'Parent',axes1);
% xlim([tps tpe2]);
% ylim([-200 650]);
% ylabel('V_{inj} [V]', 'FontSize', 16);
% legend('X-inj', 'Y-inj', 'Location', 'SouthWest', 'Orientation', 'Horizontal');
% text(tpe2 - .1, -200 + .9*(650 - -200), '(a)', 'FontSize', 14);
% text(1.51, -100, int2str(shot), 'Color', 'k', 'Fontsize', 14);
% 
% axes2 = axes('Parent',h1,'YTick',ytick_a2,'XTickLabel','',...
%     'XTick',x_ticks,...
%     'Position',[0.15 0.625 0.84 0.17],...
%     'FontSize',16);
% box('on');
% grid('on');
% hold('all');
% plot(t.psixinj*1e3, sihi_psixinj*1e3, 'r', 'LineWidth', 2,'Parent',axes2);
% plot(t.psiyinj*1e3, sihi_psiyinj*1e3, '--b', 'LineWidth', 2,'Parent',axes2);
% xlim([tps tpe2]);
% ylim(ytick_l2);
% text(tpe2 - .1, ytick_l2(1) + .9*(ytick_l2(2) - ytick_l2(1)), '(b)', 'FontSize', 14);
% ylabel('\Psi_{inj} [mWb]', 'FontSize', 16);
% 
% axes3 = axes('Parent',h1,'YTick',ytick_a3,'XTickLabel','',...
%     'XTick',x_ticks,...
%     'Position',[0.15 0.45 0.84 0.17],...
%     'FontSize',16);
% box('on');
% grid('on');
% hold('all');
% plot(t.ixinj*1e3, sihi_ixinj*1e-3, 'r', 'LineWidth', 2,'Parent',axes3);
% plot(t.iyinj*1e3, sihi_iyinj*1e-3, '--b', 'LineWidth', 2,'Parent',axes3);
% plot(t.itor*1e3, sig.itor*1e-3, 'k', 'LineWidth', 2,'Parent',axes3);
% legend('X-inj', 'Y-inj', 'I_{tor}', 'Location', 'NorthWest');
% xlim([tps tpe2]);
% ylim(ytick_l3);
% text(tpe2 - .1, ytick_l3(1) + .9*(ytick_l3(2) - ytick_l3(1)), '(c)', 'FontSize', 14);
% ylabel('I [kA]', 'FontSize', 16);
% 
% axes4 = axes('Parent',h1,'YTick',ytick_a4,'XTickLabel','',...
%     'XTick',x_ticks,...
%     'Position',[0.15 0.275 0.84 0.17],...
%     'FontSize',16);
% box('on');
% grid('on');
% hold('all');
% plot(t.den(I_den:end)*1e3, sig.den(I_den:end)*1e-19, 'k', 'LineWidth', 2,'Parent',axes4);
% xlim([tps tpe2]);
% ylim(ytick_l4);
% text(tpe2 - .1, ytick_l4(1) + .9*(ytick_l4(2) - ytick_l4(1)), '(d)', 'FontSize', 14);
% ylabel('\langlen_e\rangle [10^{19} m^{-3}]','FontSize',16);
% 
% [aa, Is] = min((tbase - s_ratio).^2);
% [aa, Ie] = min((tbase - e_ratio).^2);
% 
% axes5 = axes('Parent',h1,'YTick',ytick_a5,...
%     'XTick',x_ticks,...
%     'Position',[0.15 0.1 0.84 0.17],...
%     'FontSize',16);
% box('on');
% grid('on');
% hold('all');
% plot(tbase(Is: Ie)*1e3, i_ratio(Is: Ie), 'k', 'LineWidth', 2,'Parent',axes5);
% xlim([tps tpe2]);
% ylim(ytick_l5);
% text(tpe2 - .1, ytick_l5(1) + .9*(ytick_l5(2) - ytick_l5(1)), '(e)', 'FontSize', 14);
% ylabel('I_{tor} / I_{inj}', 'FontSize', 16);
% xlabel('t [ms]', 'FontSize', 18);
% 

%fig_save = getframe(h1);
%[Xfig, mapfig] = frame2im(fig_save);
%imwrite(Xfig, [newdir '\params_ratio' int2str(shot) date '.png']);





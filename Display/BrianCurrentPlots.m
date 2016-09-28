%% this program compares the IMP data to Taylor

clear all;
import MDSplus.*
addpath('S:\Matlab\filters');

newdir = ['S:\Matlab\high_freq_pop\figures\parameters'];
mkdir(newdir);

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

shot = 122385; % input('What shot? ');
[a, status] = mdsopen('landau.hit::hitsi', shot);

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
%% this brings in the data

[t.vxinj, dt.vxinj, sig.vxinj] = gen_data_in('\v_inj_x');
[t.vyinj, dt.vyinj, sig.vyinj] = gen_data_in('\v_inj_y');
[t.psixinj, dt.psixinj, sig.psixinj] = gen_data_in('\psi_inj_x');
[t.psiyinj, dt.psiyinj, sig.psiyinj] = gen_data_in('\psi_inj_y');
[t.ixinj, dt.ixinj, sig.ixinj] = gen_data_in('\i_inj_x');
[t.iyinj, dt.iyinj, sig.iyinj] = gen_data_in('\i_inj_y');
[t.itor, dt.itor, sig.itor] = gen_data_in('\i_tor_spaavg');
[t.den, dt.den, sig.den] = gen_data_in_smart('\n_avg_s1');

[t.itorsm, dt.itorsm, sig.itorsm] = gen_data_in('sihi_smooth(\i_tor_spaavg)');
[t_qinj, dt_qinj, sig_qinj] = ...
    gen_data_in('sihi_smooth(sigadd_quad(\i_inj_x,\i_inj_y))');
int_itor = interp1(t.itorsm, sig.itorsm, tbase);
int_iquad = interp1(t_qinj, sig_qinj, tbase);
i_ratio = int_itor./int_iquad;

% calculate the signal amplitudes here
hil_ixinj = hilbert_bsv(sig.ixinj');
hil_iyinj = hilbert_bsv(sig.iyinj');
hil_psixinj = hilbert_bsv(sig.psixinj');
hil_psiyinj = hilbert_bsv(sig.psiyinj');
hil_vxinj = hilbert_bsv(sig.vxinj');
hil_vyinj = hilbert_bsv(sig.vyinj');

mag_ixinj = sqrt(sig.ixinj.^2 + hil_ixinj'.^2);
mag_iyinj = sqrt(sig.iyinj.^2 + hil_iyinj'.^2);
mag_psixinj = sqrt(sig.psixinj.^2 + hil_psixinj'.^2);
mag_psiyinj = sqrt(sig.psiyinj.^2 + hil_psiyinj'.^2);
mag_vxinj = sqrt(sig.vxinj.^2 + hil_vxinj'.^2);
mag_vyinj = sqrt(sig.vyinj.^2 + hil_vyinj'.^2);

% smooth the signals over an injector cycle
mdsvalue('_ixinj = make_signal($1, *, $2)', mag_ixinj, t.ixinj);
mdsvalue('_iyinj = make_signal($1, *, $2)', mag_iyinj, t.iyinj);
mdsvalue('_psixinj = make_signal($1, *, $2)', mag_psixinj, t.psixinj);
mdsvalue('_psiyinj = make_signal($1, *, $2)', mag_psiyinj, t.psiyinj);
mdsvalue('_vxinj = make_signal($1, *, $2)', mag_vxinj, t.vxinj);
mdsvalue('_vyinj = make_signal($1, *, $2)', mag_vyinj, t.vyinj);

% sihi_ixinj = mdsvalue('data(sihi_smooth(_ixinj))');
% sihi_iyinj = mdsvalue('data(sihi_smooth(_iyinj))');
% sihi_psixinj = mdsvalue('data(sihi_smooth(_psixinj))');
% sihi_psiyinj = mdsvalue('data(sihi_smooth(_psiyinj))');
% sihi_vxinj = mdsvalue('data(sihi_smooth(_vxinj))');
% sihi_vyinj = mdsvalue('data(sihi_smooth(_vyinj))');

sihi_ixinj = mag_ixinj;
sihi_iyinj = mag_iyinj;
sihi_psixinj = mag_psixinj;
sihi_psiyinj = mag_psiyinj;
sihi_vxinj = mag_vxinj;
sihi_vyinj = mag_vyinj;

mdsclose;
%% plotting routines here

ts_den = -.2;
% ts_den = .00092; %start plot from here b/c of fringe jumps
[aa, I_den] = min((t.den - ts_den).^2);

if shot == 129175
    tps = -0.01;
    tpe2 = 2.3;

    ytick_a2 = [0 .25 .5];
    ytick_l2 = [-.1 .6];
    
%     ytick_a3 = [0 4 8];
%     ytick_l3 = [-1 10];
    ytick_a3 = [0 15 30];
    ytick_l3 = [-5 31];
    
    ytick_a4 = [0 .9 1.8];
    ytick_l4 = [-.2 2];
    ytick_a5 = [0 1.5 3];
    ytick_l5 = [-1 4];
    s_ratio = .000683;
    e_ratio = .00183;

elseif shot == 128436
    tps = -0.01;
    tpe2 = 2.5;

    ytick_a2 = [0 .25 .5];
    ytick_l2 = [-.1 .6];
    
%     ytick_a3 = [0 5 10];
%     ytick_l3 = [-1 12];
    ytick_a3 = [0 15 30];
    ytick_l3 = [-5 31];
    
    ytick_a4 = [0 1.5 3];
    ytick_l4 = [-.5 4];
    ytick_a5 = [0 1.5 3];
    ytick_l5 = [-1 4];
    s_ratio = .00067;
    e_ratio = .00198;
    
elseif shot == 126482
    tps = -0.01;
    tpe2 = 2.5;

%     ytick_a2 = [0 .5 1];
%     ytick_l2 = [-.1 1.2];
    ytick_a2 = [0 .5 1];
    ytick_l2 = [-.1 1.2];
    
    ytick_a3 = [0 20 40];
    ytick_l3 = [-8 43];
    
    ytick_a4 = [0 2 4];
    ytick_l4 = [-1 5];
    ytick_a5 = [0 1.5 3];
    ytick_l5 = [-.5 4];
    s_ratio = .000755;
    e_ratio = .00198;
    
elseif shot == 122385
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
    e_ratio = .001753;
    
end

x_ticks = [0.5 1 1.5 2];

scrsz = get(0,'ScreenSize');

% voltage flux and current
h1 = figure(8);     %[Left edge, Bottom edge, Width, Height]
set(h1, 'Position', [1 scrsz(4)/20 scrsz(3)/3 scrsz(4)/1.5],'Color',[1 1 1]);
set(h1, 'Name', 'params');
hold('all');
clf;
axes1 = axes('Parent',h1,'YTick',[0 250 500],'XTickLabel','',...
    'XTick',x_ticks,...
    'Position',[0.15 0.8 0.84 0.17],...
    'FontSize',16);
box('on');
grid('on');
hold('all');
plot(t.vxinj*1e3, sihi_vxinj, 'r', 'LineWidth', 2,'Parent',axes1);
plot(t.vyinj*1e3, sihi_vyinj, '--b', 'LineWidth', 2,'Parent',axes1);
xlim([tps tpe2]);
ylim([-200 650]);
text(tpe2 - .1, -200 + .9*(650 - -200), '(a)', 'FontSize', 14);
ylabel('V_{inj} [V]', 'FontSize', 16);
legend('X-inj', 'Y-inj', 'Location', 'SouthWest', 'Orientation', 'Horizontal');
text(1.51, -100, int2str(shot), 'Color', 'k', 'Fontsize', 14);

axes2 = axes('Parent',h1,'YTick', ytick_a2,'XTickLabel','',...
    'XTick',x_ticks,...
    'Position',[0.15 0.625 0.84 0.17],...
    'FontSize',16);
box('on');
grid('on');
hold('all');
plot(t.psixinj*1e3, sihi_psixinj*1e3, 'r', 'LineWidth', 2,'Parent',axes2);
plot(t.psiyinj*1e3, sihi_psiyinj*1e3, '--b', 'LineWidth', 2,'Parent',axes2);
xlim([tps tpe2]);
ylim(ytick_l2);
text(tpe2 - .1, ytick_l2(1) + .9*(ytick_l2(2) - ytick_l2(1)), '(b)', 'FontSize', 14);
ylabel('\Psi_{inj} [mWb]', 'FontSize', 16);

axes3 = axes('Parent',h1,'YTick',ytick_a3,'XTickLabel','',...
    'XTick',x_ticks,...
    'Position',[0.15 0.45 0.84 0.17],...
    'FontSize',16);
box('on');
grid('on');
hold('all');
plot(t.ixinj*1e3, sihi_ixinj*1e-3, 'r', 'LineWidth', 2,'Parent',axes3);
plot(t.iyinj*1e3, sihi_iyinj*1e-3, '--b', 'LineWidth', 2,'Parent',axes3);
xlim([tps tpe2]);
ylim(ytick_l3);
text(tpe2 - .1, ytick_l3(1) + .9*(ytick_l3(2) - ytick_l3(1)), '(c)', 'FontSize', 14);
ylabel('I_{inj} [kA]', 'FontSize', 16);

axes4 = axes('Parent',h1,'YTick',ytick_a4,'XTickLabel','',...
    'XTick',x_ticks,...
    'Position',[0.15 0.275 0.84 0.17],...
    'FontSize',16);
box('on');
grid('on');
hold('all');
plot(t.den(I_den:end)*1e3, sig.den(I_den:end)*1e-19, 'k', 'LineWidth', 2,'Parent',axes4);
xlim([tps tpe2]);
ylim(ytick_l4);
ylabel('\langlen_e\rangle [10^{19} m^{-3}]','FontSize',16);

axes5 = axes('Parent',h1,'YTick',[0 15 30],...
    'XTick',x_ticks,...
    'Position',[0.15 0.1 0.84 0.17],...
    'FontSize',16);
box('on');
grid('on');
hold('all');
plot(t.itor*1e3, sig.itor*1e-3, 'k', 'LineWidth', 2,'Parent',axes5);
xlim([tps tpe2]);
ylim([-7 34]);
ylabel('I_{tor} [kA]', 'FontSize', 16);
xlabel('t [ms]', 'FontSize', 18);

fig_save = getframe(h1);
[Xfig, mapfig] = frame2im(fig_save);
imwrite(Xfig, [newdir '\params' int2str(shot) '.png']);


% voltage flux and current
h1 = figure(9);     %[Left edge, Bottom edge, Width, Height]
set(h1, 'Position', [1 scrsz(4)/20 scrsz(3)/3 scrsz(4)/1.5],'Color',[1 1 1]);
set(h1, 'Name', 'params');
hold('all');
clf;
axes1 = axes('Parent',h1,'YTick',[0 250 500],'XTickLabel','',...
    'XTick',x_ticks,...
    'Position',[0.15 0.8 0.84 0.17],...
    'FontSize',16);
box('on');
grid('on');
hold('all');
plot(t.vxinj*1e3, sihi_vxinj, 'r', 'LineWidth', 2,'Parent',axes1);
plot(t.vyinj*1e3, sihi_vyinj, '--b', 'LineWidth', 2,'Parent',axes1);
xlim([tps tpe2]);
ylim([-200 650]);
ylabel('V_{inj} [V]', 'FontSize', 16);
legend('X-inj', 'Y-inj', 'Location', 'SouthWest', 'Orientation', 'Horizontal');
text(tpe2 - .1, -200 + .9*(650 - -200), '(a)', 'FontSize', 14);
text(1.51, -100, int2str(shot), 'Color', 'k', 'Fontsize', 14);

axes2 = axes('Parent',h1,'YTick',ytick_a2,'XTickLabel','',...
    'XTick',x_ticks,...
    'Position',[0.15 0.625 0.84 0.17],...
    'FontSize',16);
box('on');
grid('on');
hold('all');
plot(t.psixinj*1e3, sihi_psixinj*1e3, 'r', 'LineWidth', 2,'Parent',axes2);
plot(t.psiyinj*1e3, sihi_psiyinj*1e3, '--b', 'LineWidth', 2,'Parent',axes2);
xlim([tps tpe2]);
ylim(ytick_l2);
text(tpe2 - .1, ytick_l2(1) + .9*(ytick_l2(2) - ytick_l2(1)), '(b)', 'FontSize', 14);
ylabel('\Psi_{inj} [mWb]', 'FontSize', 16);

axes3 = axes('Parent',h1,'YTick',ytick_a3,'XTickLabel','',...
    'XTick',x_ticks,...
    'Position',[0.15 0.45 0.84 0.17],...
    'FontSize',16);
box('on');
grid('on');
hold('all');
plot(t.ixinj*1e3, sihi_ixinj*1e-3, 'r', 'LineWidth', 2,'Parent',axes3);
plot(t.iyinj*1e3, sihi_iyinj*1e-3, '--b', 'LineWidth', 2,'Parent',axes3);
plot(t.itor*1e3, sig.itor*1e-3, 'k', 'LineWidth', 2,'Parent',axes3);
legend('X-inj', 'Y-inj', 'I_{tor}', 'Location', 'NorthWest');
xlim([tps tpe2]);
ylim(ytick_l3);
text(tpe2 - .1, ytick_l3(1) + .9*(ytick_l3(2) - ytick_l3(1)), '(c)', 'FontSize', 14);
ylabel('I [kA]', 'FontSize', 16);

axes4 = axes('Parent',h1,'YTick',ytick_a4,'XTickLabel','',...
    'XTick',x_ticks,...
    'Position',[0.15 0.275 0.84 0.17],...
    'FontSize',16);
box('on');
grid('on');
hold('all');
plot(t.den(I_den:end)*1e3, sig.den(I_den:end)*1e-19, 'k', 'LineWidth', 2,'Parent',axes4);
xlim([tps tpe2]);
ylim(ytick_l4);
text(tpe2 - .1, ytick_l4(1) + .9*(ytick_l4(2) - ytick_l4(1)), '(d)', 'FontSize', 14);
ylabel('\langlen_e\rangle [10^{19} m^{-3}]','FontSize',16);

[aa, Is] = min((tbase - s_ratio).^2);
[aa, Ie] = min((tbase - e_ratio).^2);

axes5 = axes('Parent',h1,'YTick',ytick_a5,...
    'XTick',x_ticks,...
    'Position',[0.15 0.1 0.84 0.17],...
    'FontSize',16);
box('on');
grid('on');
hold('all');
plot(tbase(Is: Ie)*1e3, i_ratio(Is: Ie), 'k', 'LineWidth', 2,'Parent',axes5);
xlim([tps tpe2]);
ylim(ytick_l5);
text(tpe2 - .1, ytick_l5(1) + .9*(ytick_l5(2) - ytick_l5(1)), '(e)', 'FontSize', 14);
ylabel('I_{tor} / I_{inj}', 'FontSize', 16);
xlabel('t [ms]', 'FontSize', 18);


fig_save = getframe(h1);
[Xfig, mapfig] = frame2im(fig_save);
imwrite(Xfig, [newdir '\params_ratio' int2str(shot) date '.png']);





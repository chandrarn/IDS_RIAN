%% Settings
clear all; 
%close all; clc;
addpath('~/Magnetics/filters');
addpath('~/Magnetics/general');
addpath('~/Magnetics/BD_SP');
import MDSplus.*
%
% Frozen version used to make HIT-SI3 phasing comparison plots for CT talk.
%
% Special version to calculate delta B / B based on "injector BD modes" and
% "plasma BD modes" where the equilibrium BD mode is "B". Uses Amperian
% loop probes instead of gap probes.
%
% Modified replacing mode/injector comparison with additional mode
% subtraction. Trying to subtract equilibrium to show low fluctuations.
%
% Modified to reproduce BSV's BD plots
%
 
%% Parameters for all shots
%
%hfshots = [128436, 141120022, 150312013, 150312033];
hfshots = [129496, 160728012];
%lfshots = [129499, 160728024, 160818008, 160616012];
lfshots = [129499, 160728013];
 
n = 1;
% FIG 1
shot(n).title = ['HIT-SI, \Delta' '\phi = 90^{\circ}'];
shot(n).tLimWide = [0.4, 2.7];
shot(n).hfdenShot = [129449]; % MAIN SHOT HAS DENSITY
shot(n).lfdenShot = [129448];%[129596];
shot(n).lfdenTshift = -.486;
shot(n).exp = 'hitsi';
shot(n).lf = 14500;
shot(n).hf = 14500;
shot(n).lfkdShot = [];
 
n = 2;
% FIG 1
shot(n).title = ['HIT-SI3, \Delta' '\phi = 120^{\circ}'];
shot(n).tLimWide = [0.5, 2.8];
shot(n).hfdenShot = [];%150317025; % alternate shot number for density
shot(n).lfdenShot = [];[160728010]; % MAIN SHOT HAS DENSITY
shot(n).exp = 'analysis3';
shot(n).lf = 14500;
shot(n).hf = 14500;% Not looking at HF HITSI
shot(n).hfkdShot = [160728009]; % Kdot Replacement
 
n = 3;
% FIG 1
shot(n).title = ['HIT-SI3, \Delta' '\phi = 120^{\circ}'];
shot(n).tLimWide = [0.5, 2.8];
shot(n).hfdenShot = 150401017; % alternate shot number for density
shot(n).lfdenShot = 160803026; % alternate shot number for density
shot(n).exp = 'analysis3';
shot(n).lf = 14500;
shot(n).hf = 48200;
 
n = 4;
% FIG 1
shot(n).title = ['HIT-SI3, \Delta' '\phi = 0^{\circ}'];
shot(n).tLimWide = [0.3, 3.6];
shot(n).hfdenShot = 150331026; % NEGATIVE TOR CURRENT - FIND BETTER SHOT IF POSSIBLE
shot(n).lfdenShot = []; % MAIN SHOT HAS DENSITY
shot(n).exp = 'analysis3';
shot(n).lf = 14500;
shot(n).hf = 48200;
 
for n = 1:length(hfshots)
    %shot(n).exp = whatTree(hfshots(n));
    shot(n).xTick = [0.5 1 1.5 2 2.5 3 3.5];
    shot(n).currentsLim = [-95, 95];
    shot(n).pLim = [0, 16];
    shot(n).kDotLim = [0, 2.8];
    shot(n).denLim = [0, 8.4];
end
%}
 
%%
nPlots = 4; % plots per shot
sp_opt.sihi_smooth = 0; % must be '0' to get full signal for BD analysis
sp_opt.units = 'mks';
 
mu0 = 4 * pi * 1e-7;
a = 0.23; % minor radius, [m]
 
%% Set up Figure
fntsz = 16;
color = ['-k'; '-r'; '-g'; '-b'; '-c'; '-m';
     ':k'; ':r'; ':g'; ':b'; ':c'; ':m'];
S = get(0, 'ScreenSize');
 
figlabs = ['abcd'; 'efgh'; 'ijkl'; 'mnop'];
 
dx = 0.42;
 
dy = 0.2;%0.1;
y(1) = 0.045;
y(2) = y(1) + dy + 0.01;
y(3) = y(1) + 2*dy + 0.02;
y(4) = y(1) + 3*dy + 0.03;
% y(5) = y(1) + 4*dy + 0.08;
 
x = [0.1, 0.57, 0.07, 0.54]; % for big 2x2 arrangement
Y = [0.05,0.05];%[0.5, 0.5, 0, 0];
 
ax = [];
 
H.fig(1) = figure('Color', [1 1 1], 'Position', [0.1 * S(3), 0.1 * S(4), 0.38 * S(3), 0.72 * S(4)]); % currents, Fourier, Inj only, Plasma only
 %H.fig(1) = figure('Color', [1 1 1], 'Position', [0.1 * S(3), 0.1 * S(4), 0.35 * S(3), 0.36 * S(4)]);
for n = 1:length(hfshots)
    %% HF data in
    HitConn = Connection('landau.hit');
    HitConn.openTree(shot(n).exp, hfshots(n));
     
%     sihi_freq = data_in(HitConn, '\SIHI_FREQ');
     
    % Bring in Currents
    % Aaron does this:
    % hfitor = data_in(HitConn, '\i_tor_spaavg');
    hfitor = data_in(HitConn, '\i_tor_spaavg',0);
    hfitor.t = 1e3 * hfitor.t;
    hfitor.y = 1e-3 * hfitor.y;
    hfitor_sm = data_in(HitConn, '(\i_tor_spaavg)',shot(n).hf);
    hfitor_sm.t = 1e3 * hfitor_sm.t;
    hfitor_sm.y = 1e-3 * hfitor_sm.y;
     
    hfpinj = data_in(HitConn, '\P_INJ',0);
    hfpinj.t = 1e3 * hfpinj.t;
    hfpinj.y = 1e-6 * hfpinj.y; % W to MW
    hfpinj_sm = data_in(HitConn, '(\P_INJ)',shot(n).hf);
    hfpinj_sm.t = 1e3 * hfpinj_sm.t;
    hfpinj_sm.y = 1e-6 * hfpinj_sm.y; % W to MW
     
    hfkdot = data_in(HitConn, '\KDOT_INJ',0);
    hfkdot.t = 1e3 * hfkdot.t;
    hfkdot_sm = data_in(HitConn, '\KDOT_INJ',shot(n).hf);
    hfkdot_sm.t = 1e3 * hfkdot_sm.t;
     
    if isempty(shot(n).hfdenShot) % main shot has density
        hfn = data_in(HitConn, '\N_AVG_S1',0);
        hfn.t = 1e3 * hfn.t;
        hfn.y = 1e-19 * hfn.y; % take out factor of e19
%         hfn_sm = data_in(HitConn, 'sihi_smooth(\N_AVG_S1)');
%         hfn_sm.t = 1e3 * hfn_sm.t;
%         hfn_sm.y = 1e-19 * hfn_sm.y; % take out factor of e19
    else % need alternate shot
        HitConn.closeAllTree();
        HitConn.openTree(shot(n).exp, shot(n).hfdenShot);
        hfn = data_in(HitConn, '\N_AVG_S1',0);
        hfn.t = 1e3 * hfn.t;
        hfn.y = 1e-19 * hfn.y; % take out factor of e19
%         hfn_sm = data_in(HitConn, 'sihi_smooth(\N_AVG_S1)');
%         hfn_sm.t = 1e3 * hfn_sm.t;
%         hfn_sm.y = 1e-19 * hfn_sm.y; % take out factor of e19
    end
     
    HitConn.closeAllTree();
    clear HitConn % necessary for sihi_smooth to work on hit-si and hit-si3
     
    %% LF data in
    HitConn = Connection('landau.hit');
    HitConn.openTree(shot(n).exp, lfshots(n));
     
%     sihi_freq = data_in(HitConn, '\SIHI_FREQ');
     
    % Bring in Currents
    lfitor = data_in(HitConn, '\i_tor_spaavg',0);
    lfitor.t = 1e3 * lfitor.t;
    lfitor.y = 1e-3 * lfitor.y;
    lfitor_sm = data_in(HitConn, '(\i_tor_spaavg)',shot(n).lf);
    lfitor_sm.t = 1e3 * lfitor_sm.t;
    lfitor_sm.y = 1e-3 * lfitor_sm.y;
     
    lfpinj = data_in(HitConn, '\P_INJ',0);
    lfpinj.t = 1e3 * lfpinj.t;
    lfpinj.y = 1e-6 * lfpinj.y; % W to MW
    lfpinj_sm = data_in(HitConn, '(\P_INJ)',shot(n).lf);
    lfpinj_sm.t = 1e3 * lfpinj_sm.t;
    lfpinj_sm.y = 1e-6 * lfpinj_sm.y; % W to MW
     
    lfkdot = data_in(HitConn, '\KDOT_INJ',0);
    lfkdot.t = 1e3 * lfkdot.t;
    lfkdot_sm = data_in(HitConn, '(\KDOT_INJ)',shot(n).lf);
    lfkdot_sm.t = 1e3 * lfkdot_sm.t;
     
%     if n == 1 % 129499 - no density at all
%         lfn.t = [1, 2];
%         lfn.y = [5, 7]; % �\_(?)_/� -BSV
%         lfn.dy = [1, 1];
        
    if isempty(shot(n).lfdenShot) % main shot has density
        lfn = data_in(HitConn, '\N_AVG_S1',0);
        lfn.t = 1e3 * lfn.t;
        lfn.y = 1e-19 * lfn.y; % take out factor of e19
%         lfn_sm = data_in(HitConn, 'sihi_smooth(\N_AVG_S1)');
%         lfn_sm.t = 1e3 * lfn_sm.t;
%         lfn_sm.y = 1e-19 * lfn_sm.y; % take out factor of e19
    else % need alternate shot
        HitConn.closeAllTree();
        HitConn.openTree(shot(n).exp, shot(n).lfdenShot);
        lfn = data_in(HitConn, '\N_AVG_S1',0);
        lfn.t = shot(n).lfdenTshift + 1e3 * lfn.t;
        lfn.y = 1e-19 * lfn.y; % take out factor of e19
%         lfn_sm = data_in(HitConn, 'sihi_smooth(\N_AVG_S1)');
%         lfn_sm.t = 1e3 * lfn_sm.t;
%         lfn_sm.y = 1e-19 * lfn_sm.y; % take out factor of e19
    end
    
    if ~isempty(shot(n).hfkdShot) % main shot has density
         HitConn.closeAllTree();
        HitConn.openTree(shot(n).exp, shot(n).hfkdShot);
       hfkdot = data_in(HitConn, '\KDOT_INJ',0);
        hfkdot.t = 1e3 * hfkdot.t;
        hfkdot_sm = data_in(HitConn, '(\KDOT_INJ)',shot(n).hf);
        hfkdot_sm.t = 1e3 * hfkdot_sm.t;

    end
     
    HitConn.closeAllTree();
    clear HitConn % necessary for sihi_smooth to work on hit-si and hit-si3
     
%     if strcmp('hitsi3', whatTree(shots(n)))
%         iinja = data_in(HitConn, '\i_inj_a');
%         iinja.t = 1e3 * iinja.t;
%         iinja.y = 1e-3 * iinja.y;
%         iinjb = data_in(HitConn, '\i_inj_b');
%         iinjb.t = 1e3 * iinjb.t;
%         iinjb.y = 1e-3 * iinjb.y;
%         iinjc = data_in(HitConn, '\i_inj_c');
%         iinjc.t = 1e3 * iinjc.t;
%         iinjc.y = 1e-3 * iinjc.y;
%         
%         % Manually calculate quad
%         iinjq = data_in(HitConn, 'sqrt( sigadd( sigadd( sigmul(\i_inj_a, \i_inj_a), sigmul(\i_inj_b, \i_inj_b) ), sigmul(\i_inj_c, \i_inj_c) ) )');
%         iinjq.t = 1e3 * iinjq.t;
%         iinjq.y = 1e-3 * iinjq.y;
%         iinjq_sm = data_in(HitConn, 'sihi_smooth( sqrt( sigadd( sigadd( sigmul(\i_inj_a, \i_inj_a), sigmul(\i_inj_b, \i_inj_b) ), sigmul(\i_inj_c, \i_inj_c) ) ) )');
%         iinjq_sm.t = 1e3 * iinjq_sm.t;
%         iinjq_sm.y = 1e-3 * iinjq_sm.y;
%     else
%         iinjx = data_in(HitConn, '\i_inj_x');
%         iinjx.t = 1e3 * iinjx.t;
%         iinjx.y = 1e-3 * iinjx.y;
%         iinjy = data_in(HitConn, '\i_inj_y');
%         iinjy.t = 1e3 * iinjy.t;
%         iinjy.y = 1e-3 * iinjy.y;
%         
%         % Calculate Quad Currents
%         iinjy.y2 = interp1(iinjy.t, iinjy.y, iinjx.t); % put y inj current on x timebase
%         iinjq.t = iinjx.t;
%         iinjq.y = sqrt(iinjx.y.^2 + iinjy.y2.^2);
%         
%         iinjq_sm.y = sihi_smooth(iinjq.y, 1e-3 * iinjq.t, sihi_freq);
%         iinjq_sm.t = iinjq.t;
%     end      
%     end
         
%     %% ---------------------------------------------------------------- %%
     
    %% Plotting Fig 
    figure(H.fig(1)) % make current
     
    ax(length(ax) + 1) = subplot('position', [x(n), Y(n) + y(4), dx, dy]);
    set(gca, 'FontSize', fntsz);
    hold on; box on; grid on;
    H.p(1) = plot(lfitor.t, lfitor.y, '-r', 'LineWidth', 1.5);
    H.p(2) = plot(hfitor.t, hfitor.y, '-b', 'LineWidth', 1.5);
    plot(lfitor_sm.t, lfitor_sm.y, '-k', 'LineWidth', 1.5);
    plot(hfitor_sm.t, hfitor_sm.y, '-k', 'LineWidth', 1.5);
   
    legend(H.p, num2str(lfshots(n)), num2str(hfshots(n)));
    set(gca, 'XLim', shot(n).tLimWide, 'XTickLabel', [], 'XTick', shot(n).xTick);
    text(0.1, 0.84, 'I_{TOR}', 'Units', 'normalized', 'FontSize', fntsz);
    if or(n == 1, n == 3)
        ylabel('[kA]');
    else
        set(gca, 'YTickLabel', []);
    end
    set(gca, 'YLim', shot(n).currentsLim);
    text(0.02, 0.88, ['(' figlabs(n, 1) ')'], 'Units', 'normalized', 'FontSize', fntsz);
    title(shot(n).title);
     
 
    %% P_inj
    ax(length(ax) + 1) = subplot('position', [x(n), Y(n) + y(3), dx, dy]);
    set(gca, 'FontSize', fntsz);
    hold on; box on; grid on;
    plot(lfpinj.t, lfpinj.y, '-r', 'LineWidth', 1.5);
    plot(hfpinj.t, hfpinj.y, '-b', 'LineWidth', 1.5);
%     if n == 4
        tmp = plot(lfpinj_sm.t, lfpinj_sm.y, '-k', 'LineWidth', 1.5);
        plot(hfpinj_sm.t, hfpinj_sm.y, '-k', 'LineWidth', 1.5);
%         legend(tmp, 'f_{INJ} smoothed');
%     end
    set(gca, 'XLim', shot(n).tLimWide, 'XTickLabel', [], 'XTick', shot(n).xTick);
    set(gca, 'YLim', shot(n).pLim);
    text(0.1, 0.84, 'P_{INJ}', 'Units', 'normalized', 'FontSize', fntsz);
    if or(n == 1, n == 3)
        ylabel('[MW]');
    else
        set(gca, 'YTickLabel', []);
    end
    text(0.02, 0.88, ['(' figlabs(n, 2) ')'], 'Units', 'normalized', 'FontSize', fntsz);
     
 
    %% Kdot_inj
    ax(length(ax) + 1) = subplot('position', [x(n), Y(n) + y(2), dx, dy]);
    set(gca, 'FontSize', fntsz);
    hold on; box on; grid on;
    lfKD=plot(lfkdot.t, lfkdot.y, '-r', 'LineWidth', 1.5);
    hfKD=plot(hfkdot.t, hfkdot.y, '-b', 'LineWidth', 1.5);
%     if n == 4
        tmp = plot(lfkdot_sm.t, lfkdot_sm.y, '-k', 'LineWidth', 1.5);
        plot(hfkdot_sm.t, hfkdot_sm.y, '-k', 'LineWidth', 1.5);
%         legend(tmp, 'f_{INJ} smoothed');
%     end
    set(gca, 'XLim', shot(n).tLimWide, 'XTickLabel', [], 'XTick', shot(n).xTick);
    set(gca, 'YLim', shot(n).kDotLim);
    text(0.1, 0.84, 'dK/dt_{INJ}', 'Units', 'normalized', 'FontSize', fntsz);
    if or(n == 1, n == 3)
        ylabel('[Wb^2/s]');
    else
        set(gca, 'YTickLabel', []);
    end
    
    if n~=1 && ~isempty(shot(n).hfkdShot)
        legend(lfKD,num2str(shot(n).hfkdShot));
    end
    text(0.02, 0.88, ['(' figlabs(n, 3) ')'], 'Units', 'normalized', 'FontSize', fntsz);
 
    %% Density
    ax(length(ax) + 1) = subplot('position', [x(n), Y(n) + y(1), dx, dy]);
    set(gca, 'FontSize', fntsz);
    hold on; box on; grid on;
%     if n == 1
%         errorbar(lfn.t, lfn.y, lfn.dy, '-r', 'LineWidth', 1.5);
%         legend('estimated');
%     else
        lfl = plot(lfn.t, lfn.y, '-r', 'LineWidth', 1.5);
%     end
    hfl = plot(hfn.t, hfn.y, '-b', 'LineWidth', 1.5);
   % if n ~= 1
        if and(~isempty(shot(n).lfdenShot), ~isempty(shot(n).hfdenShot))
            legend([lfl, hfl], num2str(shot(n).lfdenShot), num2str(shot(n).hfdenShot));
        elseif ~isempty(shot(n).lfdenShot)
            legend(lfl, num2str(shot(n).lfdenShot));
        elseif ~isempty(shot(n).hfdenShot)
            legend(hfl, num2str(shot(n).hfdenShot));
        end
    %end
         
    set(gca, 'XLim', shot(n).tLimWide, 'XTick', shot(n).xTick);
    set(gca, 'YLim', shot(n).denLim);
    text(0.1, 0.84, 'n_e', 'Units', 'normalized', 'FontSize', fntsz);
    if or(n == 1, n == 3)
        ylabel('\times 10^{19} [m^{-3}]');
    else
        set(gca, 'YTickLabel', []);
    end
    text(0.02, 0.88, ['(' figlabs(n, 4) ')'], 'Units', 'normalized', 'FontSize', fntsz);
     
    xlabel('Time [ms]');
     
 
    %%
    %}
     
     
end
 
% set(ax, 'LineWidth', 1.5)

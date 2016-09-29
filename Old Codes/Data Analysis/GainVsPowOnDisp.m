% program to generate triplets of sigma velocity(displacement), Current gain, and
% injector power. These points will be the average of one injector cycle.
% The program will output a struct of shotnumbers, containing the data.

% load shot, generate current gain and power vectors, on exactly the same
% time base, pick same ten points from each set, average power and gain for
% the period, get sig V, store in first position. Also store time midpoint
% in another element in struct

%Note: scattered interpolant: to produce the rest of the grid
% Note: will have to squeeze the X and Y grids into column vectors, to get
% the correspoinding points from, then squeeze back into a square.
import MDSplus.*

shots = [129810,129817,129819,129820,  129495, 129496, 129793,   129440, 129441, 129443, 129446, 129449, 129450, 129451];
tLim = [1.46,1.85;1.28,1.8;1.46,1.84;1.32,1.74;  1.57,1.94; 1.51,1.94; 1.4,1.85;   1.74,1.95; 1.54,1.94; 1.34,1.94; 1.61,1.93; 1.66,1.95; 1.68,1.94; 1.41,1.95 ];
HitTree = Connection('landau.hit');

cd('T:\IDS\Data Repository\TEMP');

torPlot = 0;
chan_ranget = [8:27]; % toroidal, mohawk port in midplane
%chan_rangep = [47:61]; % poloidal
chan_rangep = [50:58]; %poloidal centered on spheromak

if torPlot
    chan_range = chan_ranget;
else
    chan_range = chan_rangep;
end

%TimeBase for CurrAmp
ts = -0.0001; % start time for new time base
%ts2 = ts - 0.0002; % start time for cropping data in TEST LATER
te = .005; % end time for new time base
%te2 = te + 0.0002; % end time for cropping data in
dtn = 0.4e-6; % dt for new time base (2.5 MHz)
tbase = ts: dtn: te; % new time base

data = zeros(length(shots)*10,5);
ind = 1;
for i = 1:length(shots)
    shots(i)
    
    cd('T:\IDS\Data Repository\TEMP');
    eval(sprintf('load(''dat%0.0f10.mat'');', shots(i))) % Load data
    cd('T:\IDS\Display');
    dat = trimRandT(dat,chan_range,tLim(i,:)); % trim to window
   
    HitTree.openTree('hitsi',shots(i));% open shot
   
    % Get Current Amplification
    [t_itorsm, dt_itorsm,itorsm] = gen_data_in('sihi_smooth(\i_tor_spaavg)',HitTree);
    [t_qinj,dt_qinj, qinj] = gen_data_in('sihi_smooth(sigadd_quad(\i_inj_x,\i_inj_y))',HitTree);
    int_itor = interp1(t_itorsm, itorsm, tbase);
    int_iquad = interp1(t_qinj, qinj, tbase);
    i_ratio = abs(int_itor./int_iquad);
    % Get Power
    [t_pinj,dt_pinj, pinj] = gen_data_in('sihi_smooth(sigadd_quad(\p_inj_x,\p_inj_y))',HitTree);
    p_inj = interp1(t_pinj,pinj,tbase);
    
    % time point
    tpoint = 1;
    
    % begin while loop
    while tpoint <=length(dat.time)-10
        sigV = std(mean(dat.vel(tpoint:tpoint+9,:)));
        
        % Gain
        [aa, Is] = min((tbase - dat.time(tpoint)*1e-3).^2);
        [aa, Ie] = min((tbase - dat.time(tpoint+9)*1e-3).^2);
        CurrAmp = mean(i_ratio(Is:Ie));
        %Power
        InjPow = mean(p_inj(Is:Ie));
        % Store in array
        data(ind,:)=[sigV,CurrAmp,InjPow*1e-6,dat.time(tpoint+5),shots(i)];
        %increment
        ind = ind+1;
        tpoint = tpoint + 10;
    end
end

%put in array
sigV = data(1:ind-1,1);
gain = data(1:ind-1,2);
powr = data(1:ind-1,3);
%make grid
[xq,yq]=meshgrid(2:.001:3,8:.0055:13.5);
vi=griddata(gain,powr,sigV,xq,yq);
figure; surf(xq,yq,vi); shading interp;

% sorted
[A,I] = sort(gain);
gainSorted = [sigV(I),gain(I)];
[B,I2] = sort (powr);
powrSorted = [sigV(I2),powr(I2)];

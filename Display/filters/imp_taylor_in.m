function [btay] = ...
    imp_taylor_in(tbase, Tscale, Sph_scale, inj_scale, ...
    mr_rad, mr_pol, mr_tor)

% this brings in all of the Taylor internal probe data

% bring in data for Taylor
[t.itor, dt.itor, sig.itor] = gen_data_in('\i_tor_spaavg');
[t.ixinj, dt.ixinj, sig.ixinj] = gen_data_in('sub_baseline_string("\\i_inj_x")');
[t.iyinj, dt.iyinj, sig.iyinj] = gen_data_in('sub_baseline_string("\\i_inj_y")');

% interpolating all data onto the same time base
int.ixinj = interp1(t.ixinj, sig.ixinj, tbase);
int.iyinj = interp1(t.iyinj, sig.iyinj, tbase);
int.itor = interp1(t.itor, sig.itor, tbase);

% ------------------------------------------------- %
% bringing in the data files for the Taylor fields

% the mid stem
% the spheromak field
s_mid = importdata('S:\Matlab\taylordata\chord.S_mid.dat', ' ', 1);
mid_x = s_mid.data(:,1); % x location
mid_y = s_mid.data(:,2); % y location
mid_z = s_mid.data(:,3); % z location
mid_r = -sqrt(2)/2*(mid_x + mid_y); % major radius
s_mid.rad = s_mid.data(:,4)*Tscale*Sph_scale;
s_mid.tor = s_mid.data(:,5)*Tscale*Sph_scale;
s_mid.pol = s_mid.data(:,6)*Tscale*Sph_scale;

% the x-inj field
x_mid = importdata('S:\Matlab\taylordata\chord.X_mid.dat', ' ', 1);
% the location data is the same as s_mid
x_mid.rad = x_mid.data(:,4)*Tscale;
x_mid.tor = x_mid.data(:,5)*Tscale;
x_mid.pol = x_mid.data(:,6)*Tscale;

% the y-inj field
y_mid = importdata('S:\Matlab\taylordata\chord.Y_mid.dat', ' ', 1);
% the location data is the same as s_mid
y_mid.rad = y_mid.data(:,4)*Tscale;
y_mid.tor = y_mid.data(:,5)*Tscale;
y_mid.pol = y_mid.data(:,6)*Tscale;

% -------------------------------------------------------- %
% scaling the Taylor equilibrium fields
% based upon the respective currents

btay_imp.rad = x_mid.rad*int.ixinj*inj_scale ...
    + y_mid.rad*int.iyinj*inj_scale ...
    + s_mid.rad*int.itor;

btay_imp.pol = x_mid.pol*int.ixinj*inj_scale ...
    + y_mid.pol*int.iyinj*inj_scale ...
    + s_mid.pol*int.itor;

btay_imp.tor = x_mid.tor*int.ixinj*inj_scale ...
    + y_mid.tor*int.iyinj*inj_scale ...
    + s_mid.tor*int.itor;

% more sign corrections
% b/c our sign conventions for positive poloidal fields are
% in the opposite direction
btay_imp.pol = -btay_imp.pol;

% find points closest to the probe points
for j = 1: length(mr_rad) % radial
    [aa, int.pr(j)] = min((mid_r - mr_rad(j)).^2);
end
for j = 1: length(mr_pol) % poloidal
    [aa, int.pp(j)] = min((mid_r - mr_pol(j)).^2);
end
for j = 1: length(mr_tor) % toroidal
    [aa, int.pt(j)] = min((mid_r - mr_tor(j)).^2);
end

% indexes of Taylor:
% probe (numbered from inner to outer most)
btay.rad = btay_imp.rad(int.pr, :);
btay.pol = btay_imp.pol(int.pp, :);
btay.tor = btay_imp.tor(int.pt, :);







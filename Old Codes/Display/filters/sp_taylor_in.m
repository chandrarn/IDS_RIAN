function [btay, phi, theta, node_save] = ...
    sp_taylor_in(tbase, sp_sph_scale, sp_ix_scale, sp_iy_scale)

% this brings in all of the Taylor surface probe data

% location of the poloidal arrays
deg = ['000'; '045'; '180'; '225'];
phi_deg = [0 45 180 225];
pol_length = 1.78568;
al_calfact = 360/pol_length;
arclength = [.177 .228 .279 .330 .622 .673 .723 .774 1.012 ...
    1.063 1.113 1.164 1.456 1.507 1.558 1.609];
th_deg = arclength*al_calfact;

% names of the probes in the poloidal arrays
rad_loc = ['S05'; 'S06'; 'S07'; 'S08'; 'L10'; 'L09'; 'L08'; 'L07'];
rad_loc_neg = ['S04'; 'S03'; 'S02'; 'S01'; 'L01'; 'L02'; 'L03'; 'L04'];

% major radius of the poloidal arrays
radius.S05 = 0.10060;     % Table of the surface probe radius for each set of probes
radius.S06 = 0.13652;     %   (this is major radius R, not r)
radius.S07 = 0.17254;     % It is necessary to look up the probe locations from the
radius.S08 = 0.20837;     % Taylor simulation files using only one dimension or there
radius.L10 = 0.43081;     % is a high probablitiy you'll overdefine the location as
radius.L09 = 0.45388;     % the machine wall is slightly different in the code than in
radius.L08 = 0.47696;     % the cad drawings.
radius.L07 = 0.50003;

% the gap probes
local_gap = ['L06'; 'L05'];
deg_gap = ['000'; '022'; '045'; '067'; '090'; '112'; '135'; '157'; ...
    '180'; '202'; '225'; '247'; '270'; '292'; '315'; '337'];
phi_deg_gap = [0 22.5 45 67.5 90 112.5 135 157.5 180 202.5 225 ...
    247.5 270 292.5 315 337.5];
arclength_gap = [.856 .930];
th_deg_gap = arclength_gap*al_calfact;

% directions of the data
dir = ['P'; 'T'];

% bring in data for Taylor
[t.itor, dt.itor, sig.itor] = gen_data_in('\i_tor_spaavg');
[t.ixinj, dt.ixinj, sig.ixinj] = gen_data_in('sub_baseline_string("\\i_inj_x")');
[t.iyinj, dt.iyinj, sig.iyinj] = gen_data_in('sub_baseline_string("\\i_inj_y")');

dead_probes = mdsvalue('\dead_probes');

% interpolating all data onto the same time base
int.ixinj = interp1(t.ixinj, sig.ixinj, tbase);
int.iyinj = interp1(t.iyinj, sig.iyinj, tbase);
int.itor = interp1(t.itor, sig.itor, tbase);

% ------------------------------------------------- %
% bringing in the data files for the Taylor fields

Taylor_Set   = 2;   % 1=2D Unifor Taylor, 2=3D Uniform Taylor, 3=3D Profiled Taylor
PREFIX = 'T:\';     % The root for all the paths (allows code to be run on other computers)

% the poloidal probe arrays
for j1 = 1: length(deg(:,1))
    
% Load the Injector Taylor Data
XinjFile = importdata([PREFIX 'A_Probe_Calibration\Codes- Analysis\Plasma Current\Set_' ...
    num2str(Taylor_Set) '\Standardized\' deg(j1, :) '_surface_xinj.dat']);
YinjFile = importdata([PREFIX 'A_Probe_Calibration\Codes- Analysis\Plasma Current\Set_' ...
    num2str(Taylor_Set) '\Standardized\' deg(j1, :) '_surface_yinj.dat']);
SphFile  = importdata([PREFIX 'A_Probe_Calibration\Codes- Analysis\Plasma Current\Set_' ...
    num2str(Taylor_Set) '\Standardized\' deg(j1, :) '_surface_sp.dat']);

S_Bdata  = SphFile.('data'); % Grab only the data from the SphFile structure  (skip the header)
X_Bdata  = XinjFile.('data'); % Grab only the data from the XinjFile structure (skip the header)
Y_Bdata  = YinjFile.('data'); % Grab only the data from the YinjFile structure (skip the header)

tdeg = ['tay' deg(j1, :)];

ltr.(tdeg) = S_Bdata(:, 1); % major rad location in machine
bps.(tdeg) = S_Bdata(:, 5); % poloidal field
bts.(tdeg) = S_Bdata(:, 6); % toroidal field

bpx.(tdeg) = X_Bdata(:, 5);
btx.(tdeg) = X_Bdata(:, 6);

bpy.(tdeg) = Y_Bdata(:, 5);
bty.(tdeg) = Y_Bdata(:, 6);

end

% the gap probes
for j1 = 1: length(local_gap(:,1))
    
% loading the gap data
XinjgapFile = importdata([PREFIX 'A_Probe_Calibration\Codes- Analysis\Plasma Current\Set_' ...
    num2str(Taylor_Set) '\Standardized\' local_gap(j1, :) '_surface_xinj.dat']);
YinjgapFile = importdata([PREFIX 'A_Probe_Calibration\Codes- Analysis\Plasma Current\Set_' ...
    num2str(Taylor_Set) '\Standardized\' local_gap(j1, :) '_surface_yinj.dat']);
SphgapFile  = importdata([PREFIX 'A_Probe_Calibration\Codes- Analysis\Plasma Current\Set_' ...
    num2str(Taylor_Set) '\Standardized\' local_gap(j1, :) '_surface_sp.dat']);

S_Bgapdata  = SphgapFile.('data'); % Grab only the data from the SphFile structure  (skip the header)
X_Bgapdata  = XinjgapFile.('data'); % Grab only the data from the XinjFile structure (skip the header)
Y_Bgapdata  = YinjgapFile.('data'); % Grab only the data from the YinjFile structure (skip the header)

tgapdeg = ['tay' local_gap(j1, :)];

thetagap.(tgapdeg) = S_Bgapdata(:, 3); % the toroidal location in degrees
bps.(tgapdeg) = S_Bgapdata(:, 5); % poloidal field
bts.(tgapdeg) = S_Bgapdata(:, 6); % toroidal field

bpx.(tgapdeg) = X_Bgapdata(:, 5);
btx.(tgapdeg) = X_Bgapdata(:, 6);

bpy.(tgapdeg) = Y_Bgapdata(:, 5);
bty.(tgapdeg) = Y_Bgapdata(:, 6);

end

jp = 0;
% constructing the Taylor field of the poloidal arrays
for j1 = 1: length(deg(:, 1))
    
    tdeg = ['tay' deg(j1, :)];

    for j2 = 1: length(rad_loc(:,1))

        % for probes at +z local
        halfpt = round(length(ltr.(tdeg))/2);
        
        [aa, Ipz] = min((radius.(rad_loc(j2,:)) - ...
            ltr.(tdeg)(1: halfpt)).^2);
        
        for j3 = 1: length(dir)
            node = ['B_' rad_loc(j2,:) dir(j3) deg(j1, :)];
            if sum(strcmp(node(3:end), dead_probes)) < 1
                jp = jp + 1;
                if dir(j3) == 'P'
                btay.(node) = bps.(tdeg)(Ipz)*int.itor * sp_sph_scale + ...
                    bpx.(tdeg)(Ipz)*int.ixinj * sp_ix_scale + ...
                    bpy.(tdeg)(Ipz)*int.iyinj * sp_iy_scale;
                elseif dir(j3) == 'T'
                btay.(node) = bts.(tdeg)(Ipz)*int.itor * sp_sph_scale + ...
                    btx.(tdeg)(Ipz)*int.ixinj * sp_ix_scale + ...
                    bty.(tdeg)(Ipz)*int.iyinj * sp_iy_scale;
                end
                phi.(node) = phi_deg(j1);
                theta.(node) = th_deg(j2);
                node_save(jp, :) = node;
            end
        end

        % for probes at -z local
        % the halfpt + 1 is a new addition, previous
        % attempts were off by one index
        [aa, Imzdum] = min((radius.(rad_loc(j2,:)) - ...
            ltr.(tdeg)(halfpt + 1: end)).^2);
        Imz = Imzdum + halfpt;

        for j3 = 1: length(dir)
            node = ['B_' rad_loc_neg(j2,:) dir(j3) deg(j1, :)];
            if sum(strcmp(node(3:end), dead_probes)) < 1
                jp = jp + 1;
                if dir(j3) == 'P'
                btay.(node) = bps.(tdeg)(Imz)*int.itor * sp_sph_scale + ...
                    bpx.(tdeg)(Imz)*int.ixinj * sp_ix_scale + ...
                    bpy.(tdeg)(Imz)*int.iyinj * sp_iy_scale;
                elseif dir(j3) == 'T'
                btay.(node) = bts.(tdeg)(Imz)*int.itor * sp_sph_scale + ...
                    btx.(tdeg)(Imz)*int.ixinj * sp_ix_scale + ...
                    bty.(tdeg)(Imz)*int.iyinj * sp_iy_scale;
                end
                phi.(node) = phi_deg(j1);
                theta.(node) = th_deg(j2 + 8);
                node_save(jp, :) = node;
            end
        end
    end
end

% SPs in the gap
for j2 = 1: length(local_gap(:,1))
    
    tgapdeg = ['tay' local_gap(j2,:)];
    
    for j1 = 1: length(phi_deg_gap)
        
        % searching for matching toroidal angles
        [aa, Igap] = min((phi_deg_gap(j1) - ...
            thetagap.(tgapdeg)).^2);

        for j3 = 1: length(dir)
            node = ['B_' local_gap(j2,:) dir(j3) deg_gap(j1, :)];
            if sum(strcmp(node(3:end), dead_probes)) < 1
                jp = jp + 1;
                if dir(j3) == 'P'
                btay.(node) = bps.(tgapdeg)(Igap)*int.itor * sp_sph_scale + ...
                    bpx.(tgapdeg)(Igap)*int.ixinj * sp_ix_scale + ...
                    bpy.(tgapdeg)(Igap)*int.iyinj * sp_iy_scale;
                elseif dir(j3) == 'T'
                btay.(node) = bts.(tgapdeg)(Igap)*int.itor * sp_sph_scale + ...
                    btx.(tgapdeg)(Igap)*int.ixinj * sp_ix_scale + ...
                    bty.(tgapdeg)(Igap)*int.iyinj * sp_iy_scale;
                end
                phi.(node) = phi_deg_gap(j1);
                theta.(node) = th_deg_gap(j2);
                node_save(jp, :) = node;
            end
        end
    end
end




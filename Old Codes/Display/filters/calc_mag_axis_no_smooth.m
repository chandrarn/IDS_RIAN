% this code is from Taylor
function [ma] = calc_mag_axis_no_smooth(tbase)

% ------------------------------------------------- %
% bringing in SP data for the currently open shot
[bsp_sp, phi_sp, r_sp, z_sp, node_sp] = sp_b_in_no_gap_rz(tbase);

dead_probes = mdsvalue('\dead_probes');
% this analysis requires all probes to exist
% will use neighboring probes to correct for this

if strcmp(dead_probes(1), 'L01P045') == 1
    bsp_sp.B_L01P045 = bsp_sp.B_L02P045;
    phi_sp.B_L01P045 = 45;
    r_sp.B_L01P045 = 43.08;
    z_sp.B_L01P045 = -22.08;
    node_sp = [node_sp; 'B_L01P045'];
else
    check_dead_probes;
end

if strcmp(dead_probes(2), 'S04P000') == 1
    bsp_sp.B_S04P000 = bsp_sp.B_S03P000;
    phi_sp.B_S04P000 = 0;
    r_sp.B_S04P000 = 10.06;
    z_sp.B_S04P000 = -13.04;
    node_sp = [node_sp; 'B_S04P000'];
else
    check_dead_probes;
end

if strcmp(dead_probes(3), 'L07P180') == 1
    bsp_sp.B_L07P180 = bsp_sp.B_L08P180;
    phi_sp.B_L07P180 = 180;
    r_sp.B_L07P180 = 50.00;
    z_sp.B_L07P180 = 8.50;
    node_sp = [node_sp; 'B_L07P180'];
else
    check_dead_probes;
end

% ------------------------------------------------- %
% putting all of the SPs into a matrix
bsp = zeros(length(node_sp(:,1)), length(tbase));
phi = zeros(length(node_sp(:,1)), 1);
r = zeros(length(node_sp(:,1)), 1);
z = zeros(length(node_sp(:,1)), 1);
rb = zeros(length(node_sp(:,1)), length(tbase));
zb = zeros(length(node_sp(:,1)), length(tbase));
k1 = 0;
for j1 = 1: length(node_sp(:,1))
    k1 = k1 + 1;
    bsp(k1, :) = bsp_sp.(node_sp(j1,:));
    phi(k1) = phi_sp.(node_sp(j1,:));
    r(k1) = r_sp.(node_sp(j1,:));
    z(k1) = z_sp.(node_sp(j1,:));
    rb(k1, :) = r(k1) * bsp(k1, :);
    zb(k1, :) = z(k1) * bsp(k1, :);
end

% ------------------------------------------------- %
% separating toroidal and poloidal signals
% separating the toroidal and poloidal sp signals
k1 = 0;
k2 = 0;
for j = 1: length(node_sp(:,1))
    if isfinite(1/strcmp(node_sp(j,6), 'P')) == 1
        k1 = k1 + 1;
        node_sp_pol(k1, :) = node_sp(j, :);
        bsp_pol(k1, :) = bsp(j, :);
        phi_pol(k1) = phi(j);
        r_pol(k1) = r(j);
        z_pol(k1) = z(j);
        rb_pol(k1, :) = rb(j, :);
        zb_pol(k1, :) = zb(j, :);
%     else
%         k2 = k2 + 1;
%         node_sp_tor(k2, :) = node_sp(j, :);
%         bsp_tor(k2, :) = bsp(j, :);
%         phi_tor(k2) = phi(j);
%         r_tor(k1) = r(j);
%         z_tor(k1) = z(j);
%         rb_tor(k1, :) = rb(j, :);
%         zb_tor(k1, :) = zb(j, :);
    end
end

% ------------------------------------------------- %
% separating signals at each toroidal angle

rb_pol000 = rb_pol(find(phi_pol == 0), :);
rb_pol045 = rb_pol(find(phi_pol == 45), :);
rb_pol180 = rb_pol(find(phi_pol == 180), :);
rb_pol225 = rb_pol(find(phi_pol == 225), :);

zb_pol000 = zb_pol(find(phi_pol == 0), :);
zb_pol045 = zb_pol(find(phi_pol == 45), :);
zb_pol180 = zb_pol(find(phi_pol == 180), :);
zb_pol225 = zb_pol(find(phi_pol == 225), :);

b_pol000 = bsp_pol(find(phi_pol == 0), :);
b_pol045 = bsp_pol(find(phi_pol == 45), :);
b_pol180 = bsp_pol(find(phi_pol == 180), :);
b_pol225 = bsp_pol(find(phi_pol == 225), :);

ma.r000 = sum(rb_pol000,1)./sum(b_pol000,1);
ma.r045 = sum(rb_pol045,1)./sum(b_pol045,1);
ma.r180 = sum(rb_pol180,1)./sum(b_pol180,1);
ma.r225 = sum(rb_pol225,1)./sum(b_pol225,1);

ma.z000 = sum(zb_pol000,1)./sum(b_pol000,1);
ma.z045 = sum(zb_pol045,1)./sum(b_pol045,1);
ma.z180 = sum(zb_pol180,1)./sum(b_pol180,1);
ma.z225 = sum(zb_pol225,1)./sum(b_pol225,1);

%------------------------------------------%
% find the mean and std of the mag axis
ma.ravg = mean([ma.r000; ma.r045; ma.r180; ma.r225]);
ma.rstd = std([ma.r000; ma.r045; ma.r180; ma.r225]);

ma.zavg = mean([ma.z000; ma.z045; ma.z180; ma.z225]);
ma.zstd = std([ma.z000; ma.z045; ma.z180; ma.z225]);


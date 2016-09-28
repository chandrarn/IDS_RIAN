function [bsp, phi, r, z, node_save] = sp_b_in_no_gap_rz_sihi(tbase)

% this brings in all of the Taylor surface probe data

% location of the poloidal arrays
deg = ['000'; '045'; '180'; '225'];
local = ['S05'; 'S06'; 'S07'; 'S08'; 'L10'; 'L09'; 'L08'; 'L07'; ...
    'L04'; 'L03'; 'L02'; 'L01'; 'S01'; 'S02'; 'S03'; 'S04'];
r_p = [10.06 13.65 17.25 20.84 43.08 45.39 47.70 50.00 ...
    50.00 47.70 45.39 43.08 20.84 17.25 13.65 10.06];
z_p = [13.04 16.66 20.23 23.82 22.08 17.56 13.03 8.50 ...
    -8.50 -13.03 -17.56 -22.08 -23.82 -20.23 -16.66 -13.04];
phi_deg = [0 45 180 225];

% directions of the data
dir = ['P'; 'T'];

dead_probes = mdsvalue('\dead_probes');

jp = 0;
% bringing in the data for the poloidal arrays
for j1 = 1: length(deg(:, 1))
    for j2 = 1: length(local(:, 1))
        for j3 = 1: length(dir)
            node = ['\B_' local(j2, :) dir(j3) deg(j1, :)];
            if sum(strcmp(node(4:end), dead_probes)) < 1
                jp = jp + 1;
                [tdum, dtdum, bdum] = gen_data_in(['sihi_smooth(' node ')']);

                bsp.(node(2: end)) = interp1(tdum, bdum, tbase);
                phi.(node(2: end)) = phi_deg(j1);
                r.(node(2: end)) = r_p(j2);
                z.(node(2: end)) = z_p(j2);
                node_save(jp, :) = node(2: end);
            end
        end
    end
end


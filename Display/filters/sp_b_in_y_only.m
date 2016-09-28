function [bsp, phi, theta, node_save] = sp_b_in_y_only(tbase)

% this brings in all of the Taylor surface probe data

% location of the poloidal arrays
deg = ['000'; '045'; '180'; '225'];
local = ['S05'; 'S06'; 'S07'; 'S08'; 'L10'; 'L09'; 'L08'; 'L07'];
phi_deg = [0 45 180 225];
pol_length = 1.78568;
al_calfact = 360/pol_length;
arclength = [.177 .228 .279 .330 .622 .673 .723 .774];
th_deg = arclength*al_calfact;

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
                [tdum, dtdum, bdum] = gen_data_in(node);

                bsp.(node(2: end)) = interp1(tdum, bdum, tbase);
                phi.(node(2: end)) = phi_deg(j1);
                theta.(node(2: end)) = th_deg(j2);
                node_save(jp, :) = node(2: end);
            end
        end
    end
end


function [bsp, phi, theta, node_save] = sp_b_in_updated(tbase)
% this brings in all of the Taylor surface probe data

% Written by Brian Victor
% Updated to newMDS by Rian Chandra

deadNans = 1; % turn dead probes into nans ( useful for symetrical plots )
              % only effects bsp struct, cause that's all I care about

% shot you want:
shotnum = 129817;

% Connect to the tree
import MDSplus.*
HitConn = Connection('landau.hit');
HitConn.openTree('analysis',shotnum);

%gen Data In
cd('T:\IDS\Display');

% location of the poloidal arrays
deg = ['000'; '045'; '180'; '225'];
local = ['S05'; 'S06'; 'S07'; 'S08'; 'L10'; 'L09'; 'L08'; 'L07'; ...
    'L04'; 'L03'; 'L02'; 'L01'; 'S01'; 'S02'; 'S03'; 'S04'];
phi_deg = [0 45 180 225];
pol_length = 1.78568;
al_calfact = 360/pol_length;
arclength = [.177 .228 .279 .330 .622 .673 .723 .774 1.012 ...
    1.063 1.113 1.164 1.456 1.507 1.558 1.609];
th_deg = arclength*al_calfact;

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

% Dead probes:
HitTree = Tree('analysis',shotnum); % This is a little silly, but the Connection object has trouble 
dead_probes = NATIVEvalue(HitTree.getNode('\dead_probes').getData())% returning this, for some reason


jp = 0;
% bringing in the data for the poloidal arrays
for j1 = 1: length(deg(:, 1))
    for j2 = 1: length(local(:, 1))
        for j3 = 1: length(dir)
            node = ['\B_' dir(j3) '_' local(j2, :) '_' deg(j1, :)];
            nodeCompare = ['\B_' local(j2, :) dir(j3) deg(j1, :)];
            if isempty(strmatch(nodeCompare(4:end), dead_probes))
                jp = jp + 1;
                [tdum, dtdum, bdum] = gen_data_in(node,HitConn);

                bsp.(node(2: end)) = interp1(tdum, bdum, tbase);
                phi.(node(2: end)) = phi_deg(j1);
                theta.(node(2: end)) = th_deg(j2);
                node_save(jp, :) = node(2: end);
            elseif deadNans % if it finds a dead node: NANS EVERYWHERE
                bsp.(node(2: end)) = ones(1,2251)*NaN;
            end
        end
    end
end

% bringing in the data for the gap arrays
for j1 = 1: length(local_gap(:, 1))
    for j2 = 1: length(deg_gap(:, 1))
        for j3 = 1: length(dir)
            node = ['\B_' dir(j3) '_' local_gap(j1, :) '_' deg_gap(j2, :)];
            nodeCompare = ['\B_' local_gap(j1, :) dir(j3) deg_gap(j2, :)];
            if isempty(strmatch(nodeCompare(4:end), dead_probes))
                jp = jp + 1;
                [tdum, dtdum, bdum] = gen_data_in(node,HitConn);
                bsp.(node(2: end)) = interp1(tdum, bdum, tbase);
                phi.(node(2: end)) = phi_deg_gap(j2);
                theta.(node(2: end)) = th_deg_gap(j1);
                node_save(jp, :) = node(2: end);
            end
        end
    end
end
end


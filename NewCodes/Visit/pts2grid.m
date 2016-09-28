function [pts, vec] = pts2grid(pts, vec, tor, config, gridShift, lengthBeyond)
%% Preamble
% This takes in points and unit vectors representing IDS chord averaged
% measurements and makes a false grid by making two copies of the points
% and vectors on either side of the original plane.
%
% 'tor' specifies which fiber sub bundle.  1 = 1:36, 2 = 37:72.
% 'config' specifies the configuration (where the lens is on the machine)
% 'dat' is included because it contains 'dat.peaks' which specifies the
% channel range EXCLUDING DEAD CHANNELS (very important)
%
% 'pts' is coordinates (n points) x (x y z)
% 'vec' is normalized vectors (n points) x (x y z)
% NB: sign convention is POSITIVE vector TOWARD the lens
%
%% Input Settings
% Everything the user might change
%
testing = 0; % plot 'pts' and 'vec' to confirm correct output

%% Calculations
% Save initial 'pts' and 'vec' for testing plotting.
%
pts_o = pts;
vec_o = vec;

vec = vertcat(vec, vec); % same for all configurations

if tor
    switch config
        case {1, 2} % in toroidal midplane, offset in +/- z
            pts2 = pts;
            pts2(:, 3) = pts2(:, 3) + gridShift;
            pts3 = pts;
            pts3(:, 3) = pts3(:, 3) - gridShift;
            pts = vertcat(pts2, pts3);
            
        case 3 % toroidal, axial port at 135 degrees
            phi = 135;
            
            dx = gridShift * sind(phi);
            dy = gridShift * cosd(phi);
            
            pts2 = pts;
            pts2(:, 1) = pts(:, 1) + dx;
            pts2(:, 2) = pts(:, 2) - dy;
            
            pts3 = pts;
            pts3(:, 1) = pts(:, 1) - dx;
            pts3(:, 2) = pts(:, 2) + dy;
            
            pts = vertcat(pts2, pts3);
            
        case 4 % toroidal, mohawk port perp.

    end
else % analyzing 'poloidal fiber' only, (maximum 37:72)
    
    switch config 
        case 1 % poloidal at toroidal angle 45 degrees
            phi = 45;
        case 2 % poloidal at toroidal angle 270 degrees
            phi = 270;
        case 3 % poloidal at toroidal angle 135 degrees
            phi = 135;
    end
    
    dx = gridShift * sind(phi);
    dy = gridShift * cosd(phi);
    
    pts2 = pts;
    pts2(:, 1) = pts(:, 1) + dx;
    pts2(:, 2) = pts(:, 2) - dy;
    
    pts3 = pts;
    pts3(:, 1) = pts(:, 1) - dx;
    pts3(:, 2) = pts(:, 2) + dy;
    
    pts = vertcat(pts2, pts3);
    
end

if testing % plot vectors and crash program
    figure(1)
    plot3(pts(:, 1), pts(:, 2), pts(:, 3), '*b');
    hold all;
    quiver3(pts(:, 1), pts(:, 2), pts(:, 3), ...
        lengthBeyond*vec(:, 1), lengthBeyond*vec(:, 2), lengthBeyond*vec(:, 3), ...
        'AutoScale', 'off')
    plot3(pts_o(:, 1), pts_o(:, 2), pts_o(:, 3), '*r');
    quiver3(pts_o(:, 1), pts_o(:, 2), pts_o(:, 3), ...
        lengthBeyond*vec_o(:, 1), lengthBeyond*vec_o(:, 2), lengthBeyond*vec_o(:, 3), ...
        'AutoScale', 'off')
    grid on;
    xyLim = 0.6 + lengthBeyond;
    set(gca, 'XLim', [-xyLim, xyLim], 'YLim', [-xyLim, xyLim]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    
    % do not return 'pts' and 'vec'
    disp('Crashing on purpose because in testing mode');
    clear pts vec;
end

end
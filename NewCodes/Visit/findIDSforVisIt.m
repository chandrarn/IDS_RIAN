function [pts, vec, origin] = findIDSforVisIt(tor, config, dat, ptShift, lengthBeyond)
%% Preamble
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
% Used for all configurations
%
angles = 0:2.95:2.95*35; % angles of chords
angles = angles - mean(angles); % center chords around zero angle
Rnim = 0.53385; % Machine radius [m]

%% Toroidal Fiber (channels 1:36)
%
if tor 
    switch config
        case 1 % toroidal, 71 degree port
        case 2 % toroidal, mohawk port in midplane
            
            % Locate lens first
            thetaLens = 159.4; % toroidal angle of lens in mohawk port (MACHINE coordinates)
            Rlens = .584; % radial location of lens [m]
            [Xlens, Ylens] = pol2cart(deg2rad(thetaLens), Rlens);
            origin = [Xlens, Ylens, 0]; % returned for 'makeLineOutGeom' code
            
            % Angle of Lens Axis in Machine coordinates
            aPort = 31.866; % angle of lens axis w.r.t. the radial direction!
            aRhat = thetaLens - 180; % angle of -radial direction (~ -20 deg)
            aLensAxis = aRhat + aPort; % angle of lens axis (~ +11 deg)
            
            % Angles of chords in Machine coordinates
            angles = angles + aLensAxis; % offset angles by angle of lens axis
            angles = angles(dat.peaks); % trim to only use chords in channel range
            
            % Make arrays of discrete coordinates along chords
            xarray = linspace(Xlens, Rnim + lengthBeyond, 10000);
            xpts = zeros(length(angles), 1); % preallocate
            ypts = xpts; % preallocate
            for n = 1:length(angles);
                yarray = Ylens + tand(angles(n)) * (xarray - min(xarray));
                % convert (x,y) of all chords into polar
                [~, rh] = cart2pol(xarray, yarray);
                % Find intersection of chords with boundary
                [~, ~, ~, imin] = extrema(abs(rh - Rnim));
                imin = max(imin); % keep only far side intersection index
                xmin = xarray(imin); % find 'xarray' value at intersection
                % min of [x + dx(ie: length beyond) - xarray] to find index
                [~, imin] = min(abs(xmin + lengthBeyond * cosd(angles(n)) - xarray));
                xpts(n) = xarray(imin);
                ypts(n) = yarray(imin) + ptShift;
            end
            
%             eps = 1e-4;
            eps = 0; % if zero, VisIt assumes this is 2D data
            pts = horzcat(xpts, ypts, eps*ones(length(angles), 1));
            
            % Create normalized vectors TOWARD THE LENS
            xvec = -cosd(angles);
            yvec = -sind(angles);
            eps = 0; % epsilon, so z component is not exactly zero
            zvec = eps * ones(length(angles), 1);
            vec = horzcat(xvec', yvec', zvec);
            
            if testing
                % Machine boundary circle array
                thetaNim = linspace(0, 2*pi, 10000); % array of angles in circle
                RnimVect = Rnim * ones(1, length(thetaNim));
                [x, y] = pol2cart(thetaNim, RnimVect);
                
                figure(1)
                plot(Xlens, Ylens, 'r+');
                hold on;
                plot(xpts, ypts, 'b.');
                quiver(xpts, ypts, lengthBeyond*xvec', lengthBeyond*yvec', 'AutoScale', 'off');
                plot(x, y, 'g');
                set(gca, 'XLim', [-.6 .7], 'YLim', [-.6 .7]);
            end
    
        case 3 % toroidal, axial port at 135 degrees
            
            [TH, R, Z, vecRZ, angles, rLens, zLens] = poloidalPrep(angles, dat, lengthBeyond, tor, testing);
            
            pts = zeros(length(angles), 3); % preallocate
            vec = pts; % preallocate
            
            phi = 135;
            zSign = 1;
            
            % Calculate origin (lens coordinates)
            [origin(1), origin(2), origin(3)] = pol2cart(deg2rad(phi), rLens, zLens);
            origin(3) = zSign * origin(3);
            
            TH = TH + deg2rad(phi); % shift
            [pts(:, 1), pts(:, 2), pts(:, 3)] = pol2cart(TH, R, Z);
            
            vec(:, 1) = cosd(phi) * vecRZ(:, 1); % x component
            vec(:, 2) = sind(phi) * vecRZ(:, 1); % y component
            vec(:, 3) = zSign * vecRZ(:, 2); % z mirrored from "poloidal" fiber
            
        case 4 % toroidal, mohawk port perp.

    end
else % analyzing 'poloidal fiber' only, (maximum 37:72)
    
%% Poloidal Fiber (channels 37:72)
    
    [TH, R, Z, vecRZ, angles, rLens, zLens] = poloidalPrep(angles, dat, lengthBeyond, tor, testing);
    
    pts = zeros(length(angles), 3); % preallocate
    vec = pts; % preallocate
    
    switch config % translate ptsRZ to different toroidal angles
        case 1 % poloidal at toroidal angle 45 degrees
            phi = 45;
            zSign = 1;

        case 2 % poloidal at toroidal angle 270 degrees
            phi = 270;
            zSign = 1;
            
        case 3 % poloidal at toroidal angle 135 degrees
            phi = 135;
            zSign = -1;
            
    end
    
    % Calculate origin (lens coordinates)
    [origin(1), origin(2), origin(3)] = pol2cart(deg2rad(phi), rLens, zLens);
    origin(3) = zSign * origin(3);

    TH = TH + deg2rad(phi); % shift
    Z = zSign * Z; % account for which side of machine the lens is on
    [pts(:, 1), pts(:, 2), pts(:, 3)] = pol2cart(TH, R, Z);
    
    vec(:, 1) = cosd(phi) * vecRZ(:, 1); % x component
    vec(:, 2) = sind(phi) * vecRZ(:, 1); % y component
    vec(:, 3) = zSign * vecRZ(:, 2); % z unchanged
    
end

if testing
    figure(2)
    plot3(pts(:, 1), pts(:, 2), pts(:, 3), '*b');
    hold on;
    quiver3(pts(:, 1), pts(:, 2), pts(:, 3), ...
        lengthBeyond*vec(:, 1), lengthBeyond*vec(:, 2), lengthBeyond*vec(:, 3), ...
        'AutoScale', 'off')
    plot3(origin(1), origin(2), origin(3), '*m');
    grid on;
    
    % Machine boundary circle array
    thetaNim = linspace(0, 2*pi, 10000); % array of angles in circle
    RnimVect = Rnim * ones(1, length(thetaNim));
    [xRing, yRing] = pol2cart(thetaNim, RnimVect);
    zRing = zeros(1, length(thetaNim));
    plot3(xRing, yRing, zRing, '-k');
    
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    
    % make it crash to avoid saving data
    clear vec pts
    disp('Crashed on purpose because in testing mode');
end

end

function[TH, R, Z, vecRZ, angles, rLens, zLens] = poloidalPrep(angles, dat, lengthBeyond, tor, testing)
% --- USED FOR ALL POLOIDAL PLANE CONFIGURATIONS
    
    % First find points at 0 degrees, then shift to respective angles
    rLens = 0.368; % [m], major radius of lens location
    zLens = -0.3126 - 0.0254 - 0.0127; % [m], 
    % distance from magnetic axis to wall + thickness of window + thickness of wall
    
    % Load and prepare NIMROD boundary
    load('nim_hit_coord.mat'); % contains NIMROD 'R' and 'Z' grid
    half = ceil(size(R, 2) / 2); % half the number of points along the boundary
    Rnim = [R(:, end)', R(end, half:end)];
    Znim = [Z(:, end)', Z(end, half:end)];
    
    % Make IDS chords
    zEnd = max(Znim) + lengthBeyond; % maximum extent in Z where a point could be
    zChord = linspace(zLens, zEnd, length(Rnim)); % array of 'Z' points
    zChord = -zLens + zChord; % offset so angles work
    if tor
        angles = angles(dat.peaks); % select only relevant angles
    else
        angles = angles(dat.peaks - 36); % select only relevant angles
    end
    
    rChord = zeros(length(angles), length(zChord)); % preallocate
    for n = 1:length(angles)
        rChord(n, :) = tand(angles(n)) * zChord;
    end
    rChord = rChord + rLens; % push back out to real 'R'
    zChord = zLens + zChord; % push back out to real 'Z'
    
    % Find intersection of chords with walls
    intsec = zeros(length(angles), 2); % preallocate
    ptsRZ = intsec; % preallocate
    for n = 1:length(angles)
        minVec = zeros(1, length(Rnim)); % preallocate
        for m = 1:length(Rnim)
            minVec(m) = min(abs(Rnim - rChord(n, m)) + abs(Znim - zChord(m)));
        end
        [~, imin] = min(minVec);
        intsec(n, :) = [rChord(n, imin), zChord(imin)]; % find intersection
        
        % Find points a set chord length beyond
        [~, imin] = min(abs(zChord(imin) + lengthBeyond * cosd(angles(n)) - zChord));
        ptsRZ(n, :) = [rChord(n, imin), zChord(imin)];
    end
    
    % Now find vecRZ
    vecRZ = ptsRZ; % preallocate
    vecRZ(:, 1) = -sind(angles);
    vecRZ(:, 2) = -cosd(angles);

    if testing
        figure(1)
        plot(Rnim, Znim, '*b')
        hold on;
        for n = 1:length(angles)
            plot(rChord(n, :), zChord, '*k')
        end
        plot(intsec(:, 1), intsec(:, 2), '*r')
        plot(ptsRZ(:, 1), ptsRZ(:, 2), '*r')
        quiver(ptsRZ(:, 1), ptsRZ(:, 2), lengthBeyond*vecRZ(:, 1), lengthBeyond*vecRZ(:, 2), ...
            'AutoScale', 'off')
        plot(rLens, zLens, '+r');
    end
    
    [TH, R, Z] = cart2pol(ptsRZ(:, 1), zeros(size(ptsRZ, 1), 1), ptsRZ(:, 2));
    
end
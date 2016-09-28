function [X, Y, Vx, Vy] = pol2cartVect(THETA, RHO, Vt, Vr)
%pol2cartVect 
% NB: in Matlab, the origin of THETA is at y=0, x positive.

VXr = Vr .* cos(THETA);
VYr = Vr .* sin(THETA);

VXt = -Vt .* sin(THETA);
VYt = Vt .* cos(THETA);

Vx = VXr + VXt;
Vy = VYr + VYt;

[X, Y] = pol2cart(THETA, RHO);
end


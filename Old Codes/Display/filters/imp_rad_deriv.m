function [dbrdi, dbpdi, dbtdi] = ...
    imp_rad_deriv(b, N, shot)

% this takes the derivative using data from the IMP
% the beginning and end probes use a first order deriv
% and the middle probes use a second order deriv

% b is the signal of interest
% N is the number of probes

if shot <= 118389
    % the radial distance between probes is 0.5 inches
    dr = 0.0127;
elseif shot > 118389
    % the radial distance between probes is 1 inch
    dr = 0.0254;
end

%----------------------------------------------------------------%
% b_rad derivative
dbrdi = zeros(N, length(b));

% 1st location
dbrdi(1, :) = (b(1, 2, :) - b(1, 1, :))/dr;

% 2nd order algorithm for 1st location
% dbrdi(1, :) = (-3*b(1,1,:) + 4*b(1,2,:) - b(1,3,:))/(2*dr);

% middle locations
for k = 1: N - 2
    kk = k + 1;
    dbrdi(kk, :) = (b(1, kk + 1, :) - b(1, k, :))/(2*dr);
end

% last location
dbrdi(N, :) = (b(1, N, :) - b(1, N - 1, :))/dr;

% 2nd order algorithm for last location
% dbrdi(N, :) = (b(1,N-2,:) - 4*b(1,N-1,:) + 3*b(1,N,:))/(2*dr);

%----------------------------------------------------------------%
% b_pol derivative
dbpdi = zeros(N, length(b));

% 1st location
dbpdi(1, :) = (b(2, 2, :) - b(2, 1, :))/dr;

% 2nd order algorithm for 1st location
% dbpdi(1, :) = (-3*b(2,1,:) + 4*b(2,2,:) - b(2,3,:))/(2*dr);

% middle locations
for k = 1: N - 2
    kk = k + 1;
    dbpdi(kk, :) = (b(2, kk + 1, :) - b(2, k, :))/(2*dr);
end

% last location
dbpdi(N, :) = (b(2, N, :) - b(2, N - 1, :))/dr;

% 2nd order algorithm for last location
% dbpdi(N, :) = (b(2,N-2,:) - 4*b(2,N-1,:) + 3*b(2,N,:))/(2*dr);

%----------------------------------------------------------------%
% b_tor derivative
dbtdi = zeros(N, length(b));

% 1st location
dbtdi(1, :) = (b(3, 2, :) - b(3, 1, :))/dr;

% 2nd order algorithm for 1st location
% dbtdi(1, :) = (-3*b(3,1,:) + 4*b(3,2,:) - b(3,3,:))/(2*dr);

% middle locations
for k = 1: N - 2
    kk = k + 1;
    dbtdi(kk, :) = (b(3, kk + 1, :) - b(3, k, :))/(2*dr);
end

% last location
dbtdi(N, :) = (b(3, N, :) - b(3, N - 1, :))/dr;

% 2nd order algorithm for last location
% dbtdi(N, :) = (b(3,N-2,:) - 4*b(3,N-1,:) + 3*b(3,N,:))/(2*dr);




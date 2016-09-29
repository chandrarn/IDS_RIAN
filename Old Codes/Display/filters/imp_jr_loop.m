function jr = imp_jr_loop(top, mid, bot)

% this calculates j_r using a loop and the tor and pol probes

mu0 = 4*pi*1e-7;

% l1 = l3 = 1.237 in
l1 = .0314; % m
l3 = .0314; % m

% l2 = 1.75 in
l2 = .0445; % m

% A = (1.75 in)(.875 in)/2
A = 4.94e-4; % m^2

% side 1
b1 = sqrt(2)/4*squeeze((mid(2, :, :) - mid(3, :, :) ...
    + top(2, :, :) - top(3, :, :)));

% side 2
b2 = squeeze(top(3, :, :) + bot(3, :, :))/2;

% side 3
b3 = sqrt(2)/4*squeeze((- bot(2, :, :) - bot(3, :, :) ...
    - mid(2, :, :) - mid(3, :, :)));

% int(b dl) = mu0 A j
jr = (b1*l1 + b2*l2 + b3*l3)/(mu0*A);



function jt = imp_jt_loop(top, mid, bot, N, dir, shot)

% this calculates j_t using a loop and the rad and pol probes

mu0 = 4*pi*1e-7;

if shot <= 118389
    % l1 = l3 = 0.5 in
    l1 = .0127; % m
    l3 = .0127; % m
elseif shot > 118389
    % l1 = l3 = 1 in
    l1 = .0254; % m
    l3 = .0254; % m
end

% l2 = l4 = 0.875 in
l2 = .0222; % m
l4 = .0222; % m

% A = (.5 in)(.875 in)
A = l1*l2; % m^2

for j = 1: N - 1
    
    % side 1
    b1 = squeeze(mid(1, j, :) + mid(1, j + 1, :))/2;

    % side 2
    if(strcmp(dir, 't'))
        b2 = squeeze(mid(2, j + 1, :) + top(2, j + 1, :))/2;
    elseif(strcmp(dir, 'b'))
        b2 = squeeze(mid(2, j + 1, :) + bot(2, j + 1, :))/2;
    end

    % side 3
    if(strcmp(dir, 't'))
        b3 = squeeze(- top(1, j + 1, :) - top(1, j, :))/2;
    elseif(strcmp(dir, 'b'))
        b3 = squeeze(- bot(1, j + 1, :) - bot(1, j, :))/2;
    end

    % side 4
    if(strcmp(dir, 't'))
        b4 = squeeze(- mid(2, j, :) - top(2, j, :))/2;
    elseif(strcmp(dir, 'b'))
        b4 = squeeze(- mid(2, j, :) - bot(2, j, :))/2;
    end

    % int(b dl) = mu0 A j
    jt(j, :) = (b1*l1 + b2*l2 + b3*l3 + b4*l4)/(mu0*A);
    
end




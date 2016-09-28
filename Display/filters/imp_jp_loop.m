function jp = imp_jp_loop(top, mid, bot, N, dir, shot)

% this calculates j_p using a loop and the tor and rad probes

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
    

if(strcmp(dir, 'tm'))
    % l2 = l4 = 0.875 in
    l2 = .0222; % m
    l4 = .0222; % m
elseif(strcmp(dir, 'mb'))
    % l2 = l4 = 0.875 in
    l2 = .0222; % m
    l4 = .0222; % m
elseif(strcmp(dir, 'tb'))
    % l2 = l4 = 1.75 in
    l2 = .0445; % m
    l4 = .0445; % m
end

% area of the face
A = l1*l2; % m^2

for j = 1: N - 1
    
    % side 1
    if(strcmp(dir, 'tm'))
        b1 = squeeze(- top(1, j, :) - top(1, j + 1, :))/2;
    elseif(strcmp(dir, 'mb'))
        b1 = squeeze(- mid(1, j, :) - mid(1, j + 1, :))/2;
    elseif(strcmp(dir, 'tb'))
        b1 = squeeze(- top(1, j, :) - top(1, j + 1, :))/2;
    end

    % side 2
    if(strcmp(dir, 'tm'))
        b2 = squeeze(- top(3, j + 1, :) - mid(3, j + 1, :))/2;
    elseif(strcmp(dir, 'mb'))
        b2 = squeeze(- mid(3, j + 1, :) - bot(3, j + 1, :))/2;
    elseif(strcmp(dir, 'tb'))
        b2 = squeeze(- top(3, j + 1, :) - bot(3, j + 1, :))/2;
    end

    % side 3
    if(strcmp(dir, 'tm'))
        b3 = squeeze(mid(1, j, :) + mid(1, j + 1, :))/2;
    elseif(strcmp(dir, 'mb'))
        b3 = squeeze(bot(1, j, :) + bot(1, j + 1, :))/2;
    elseif(strcmp(dir, 'tb'))
        b3 = squeeze(bot(1, j, :) + bot(1, j + 1, :))/2;
    end

    % side 4
    if(strcmp(dir, 'tm'))
        b4 = squeeze(top(3, j, :) + mid(3, j, :))/2;
    elseif(strcmp(dir, 'mb'))
        b4 = squeeze(mid(3, j, :) + bot(3, j, :))/2;
    elseif(strcmp(dir, 'tb'))
        b4 = squeeze(top(3, j, :) + bot(3, j, :))/2;
    end

    % int(b dl) = mu0 A j
    jp(j, :) = (b1*l1 + b2*l2 + b3*l3 + b4*l4)/(mu0*A);
    
end




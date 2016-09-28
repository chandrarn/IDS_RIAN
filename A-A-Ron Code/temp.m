
%% Calculate PIX_SP

% max_evals = 200;
PIX_SP = zeros(1, n_chan); % allocate array
nm = motorSpeed * 1e-3 * time; % convert 'time' to [s], [s] to [nm]
% dnm = nm(end) - nm(1); % lambda difference over whole movie

S = get(0,'ScreenSize');
fntsz = 20;
h1 = figure('Visible','on','Name','Linear Fit','Position',...
    [S(3)/10, S(4)/8, 4*S(3)/5, 3*S(4)/4], 'Color', [1 1 1]);
% store the standard deviation in fitting error
err = zeros(n_chan,1);
p_fit = zeros(n_chan,2);
Chi_sp = zeros(n_chan,1);
dPar = zeros(n_chan,2);
stddev = zeros(n_chan,size(time));
corr = zeros(n_chan,2,2);
R = zeros(n_chan,1);
for n = 1:n_chan
    % slope estimate nm/s * s/pix = [nm/pix]
%     A(1) = dnm/(centers(n, end) - centers(n, 1));
%     A(2) = nm(1); % offset estimate
%     [estimate] = optimize(centers(n, :)', nm, A, max_evals);
    centers(n, :)';
    nm;
    p = polyfit(centers(n, :), nm, 1);
    PIX_SP(n) = abs(p(1));
    err(n) = std(centers(n,:)-polyval(p,1:length(centers(n,:))));
    clf; % Plot Raw Data and Fit ----------------------------
    h4 = axes('Parent', h1, 'Position', [.2 .2 .6 .6], 'FontSize', fntsz);
    h5 = plot(centers(n, :), p(1) * centers(n, :) + p(2), '-*r');
    hold on;
    h6 = plot(centers(n, :), nm, 'ob');
    title(['Raw Wavelength vs. Pixel and Linear Fit, Channel ' num2str(PEAKS(n,1))]);
    ylabel('Wavelength [nm]');
    xlabel('Pixel (wavelength space)');
    
    % LSQ curvefitting
    guess = p;
    x = centers(n,:)';
    z = nm;
    lb = p-.001;
    ub = p+.001;
    dp = [.001,.001];
    n
     %[p_fit(n,:),Chi_sq(n),dPar(n,:),stddev(n,:),corr(n,:,:),R(n)] = ...
     %               lm(f2, guess, x, z, .001,dp);%, p_min,p_max,0)
    pause(0.1);
end
assignin('base','err',err);
assignin('base','p_fit',p_fit);
assignin('base','Chi_sq',Chi_sq);
assignin('base','dPar',dPar);
assignin('base','stddev',stddev);
assignin('base','corr',corr);
assignin('base','R',R);
err
PIX_SP = 1e-9 * PIX_SP; % convert from nm to meters

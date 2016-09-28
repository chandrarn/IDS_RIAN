function PIX_SP = fitMotor(shot3, shot4, PEAKS, motorSpeed, brightWing, xWing)
% addpath('S:\MC_IDS\Matlab Code\Core Codes\lsqcurvefit');

plotFits = 0;

% Load Data

% cd('T:\PhantomMovies');
data3 = importdata(['shot' int2str(shot3) '.mat']); % [counts] (time x wavelength space x channel space)
time3 = importdata(['t' int2str(shot3) '.mat']); % [ms]

data3 = data3(:, end:-1:1, end:-1:1);

if shot4 ~= 0 % using both fibers
    data4 = importdata(['shot' int2str(shot4) '.mat']); % [counts] (time x wavelength space x channel space)
    time4 = importdata(['t' int2str(shot4) '.mat']); % [ms]

    data4 = data4(:, end:-1:1, end:-1:1);

    % Add movies together

    if length(time3) > length(time4)
        data = data4 + data3(1:length(time4), :, :);
        time = time4;
        clear data3 data4 time3 time4;
    else
        data = data3 + data4(1:length(time3), :, :);
        time = time3;
        clear data3 data4 time3 time4;
    end
else % only using one fiber
    data = data3;
    time = time3;
    clear data3 time3;
end

% cd('S:\MC_IDS\Matlab Code\Calibration\Calibration v2');

data = cast(data, 'double');
[n_time, n_pix, n_chan] = size(data);

%% Loop over Channels

try
    options = optimsetv61('lsqcurvefit'); % set options for curve fitting
catch
    options = optimset('lsqcurvefit'); % set options for curve fitting
end
options.TolFun = 1e-8 * max(data(:)); % set tolerance for curve fitting
f = @singletGauss2D; % function handle

centers = NaN*zeros(size(PEAKS, 1), length(time)); % initialize array to store centers for all time

if plotFits
    S = get(0,'ScreenSize');
    fntsz = 20;
    h1 = figure('Visible','on','Name','Individual Fit','Position',...
        [S(3)/10, S(4)/8, 4*S(3)/5, 3*S(4)/4], 'Color', [1 1 1]);
end

for n = 1:size(PEAKS, 1)

    % find real space bound for fitting domain

    xBound = round(PEAKS(n, 2)) - xWing : round(PEAKS(n, 2)) + xWing;

    for m = 1:n_time
         if (mod(m, 10) == 0) 
            display(['CHANNEL: ' num2str(n) '/' num2str(size(PEAKS, 1)) ', TIME PT: ' num2str(m) '/' num2str(n_time)]);
         end
        % find wavelength coordinate of starting peak
        
        [amp, y] = max(max(data(m, :, xBound))); % which column has the max
        [amp, y] = max(squeeze(data(m, :, xBound(y)))); % which row in that column

        % find wavelength direction bound for fitting

        yBound = y - brightWing : y + brightWing;

        % catch in case out of bounds

        if yBound(1) < 1
            yBound = yBound(find(yBound == 1)) : yBound(end);
        end
        if yBound(end) > n_pix
            yBound = yBound(1) : yBound(find(yBound == n_pix));
        end

        % make grid

        [X, Y] = meshgrid(xBound, yBound);
        Z = squeeze(data(m, yBound, xBound));
        
        if plotFits
            clf; % Plot Raw Data ----------------------------
            ax(1) = axes('Parent', h1, 'Position', [.1 .1 .2 .8], 'FontSize', fntsz);
            h3 = surf(X, Y, Z);
            hold on;
            shading interp;
            colormap jet;
            colorbar;
            grid on;
            view([0 90]);
            set(gca, 'FontSize', fntsz);
%             set(gca, 'XLim', [xBound(1) xBound(end)], 'YLim', [yBound(1) yBound(end)]);

            xlabel('Real Space');
            ylabel('Wavelength Space');
            title(['Channel ' num2str(PEAKS(n, 1)) ', Time pt. ' num2str(m)]);
        end

        % Initial Guesses

        sigx = PEAKS(n, 4);
        sigy = PEAKS(n, 5);
        vol = 6 * sigx * sigy * amp;
        
        guess = [vol, PEAKS(n, 2), y, sigx, sigy, min(Z(:))];

        lb = [0.4*vol, guess(2)-0.5, guess(3)-0.5, 0.8*guess(4), 0.8*guess(5), 0]; % lower bound
        ub = [1.8*vol, guess(2)+0.5, guess(3)+0.5, 1.2*guess(4), 1.2*guess(5), 0.5*amp]; % upper bound

        % Reshape data and grid

        x(:, 1) = X(:);
        x(:, 2) = Y(:);
        z = Z(:);

        % Fit Gaussian ----------------------------------------
        
        parameters = lsqcurvefit(f, guess, x, z, lb, ub, options);

        centers(n, m) = parameters(3); % save center in wavelength space only

        if plotFits
            % create fine mesh --------------------------------
            nfx = 100;
            nfy = 200;
            xBoundf = linspace(xBound(1), xBound(end), nfx);
            yBoundf = linspace(yBound(1), yBound(end), nfy);
            [Xf, Yf] = meshgrid(xBoundf, yBoundf);

            % Reshape fine mesh for execution by Gaussian function

            xf(:, 1) = Xf(:);
            xf(:, 2) = Yf(:);

            zf = singletGauss2D(parameters, xf); % calculate fit on fine mesh

            Zf = reshape(zf, size(Xf, 1), size(Xf, 2)); % reshape into 2D image

            % Plot Fit ----------------------

            ax(2) = axes('Parent', h1, 'Position', [.4 .1 .2 .8], 'FontSize', fntsz);
            h4 = surf(Xf, Yf, Zf);
            hold on;
            shading interp;
            colormap jet;
            colorbar;
            grid on;
            view([0 90]);
            set(gca, 'FontSize', fntsz);
            set(gca, 'XLim', [xBound(1) xBound(end)], 'YLim', [yBound(1) yBound(end)]);

            xlabel('Real Space');
            title(['Fit, Channel ' num2str(PEAKS(n, 1)) ', Time pt. ' num2str(m)]);

            % Link Z Axis

            zmin = min(Z(:));
            zmax = max(Z(:));
            zfmin = min(Zf(:));
            zfmax = max(Zf(:));
            set(ax(1), 'ZLim', [min([zmin zfmin]), max([zmax zfmax])]);
            set(ax(2), 'ZLim', [min([zmin zfmin]), max([zmax zfmax])]);
            
            pause(0.1);

        end % plotting
        
        clear x z xf zf; % when they change size, the dimensions mismatch
        
    end % time loop
    assignin('base','centers',centers);
end % channel loop

%% Calculate PIX_SP

% max_evals = 200;
PIX_SP = zeros(1, size(PEAKS, 1)); % allocate array
nm = motorSpeed * 1e-3 * time; % convert 'time' to [s], [s] to [nm]
% dnm = nm(end) - nm(1); % lambda difference over whole movie

S = get(0,'ScreenSize');
fntsz = 20;
h1 = figure('Visible','on','Name','Linear Fit','Position',...
    [S(3)/10, S(4)/8, 4*S(3)/5, 3*S(4)/4], 'Color', [1 1 1]);

for n = 1:size(PEAKS, 1)
    % slope estimate nm/s * s/pix = [nm/pix]
%     A(1) = dnm/(centers(n, end) - centers(n, 1));
%     A(2) = nm(1); % offset estimate
%     [estimate] = optimize(centers(n, :)', nm, A, max_evals);
    p = polyfit(centers(n, :)', nm, 1);
    PIX_SP(n) = abs(p(1));
    if plotFits
        clf; % Plot Raw Data and Fit ----------------------------
        h4 = axes('Parent', h1, 'Position', [.2 .2 .6 .6], 'FontSize', fntsz);
        h5 = plot(centers(n, :), p(1)*centers(n, :) + p(2), '-r');
        hold on;
        h6 = plot(centers(n, :), nm, 'ob');
        title(['Raw Wavelength vs. Pixel and Linear Fit, Channel ' num2str(PEAKS(n,1))]);
        ylabel('Wavelength [nm]');
        xlabel('Pixel (wavelength space)');
        
        pause(0.5);
    end
end

PIX_SP = 1e-9 * PIX_SP; % convert from nm to meters

end

%% Linear fit for channel spacing
% 
% function[estimates] = optimize(pix, nm, A, max_evals)
% options = optimset('MaxFunEvals', max_evals);
% model = @linear_func;
% estimates = fminsearch(model, A);
% 
%     function [sse, nm_fit] = linear_func(A)
%         nm_fit = A(1) .* pix + A(2);
%         sse = sum((nm_fit - nm).^2);
%     end
% end








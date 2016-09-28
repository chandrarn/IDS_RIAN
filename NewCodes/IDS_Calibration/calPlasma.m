function PIX_SP = calPlasma(shotPlas, PIX_SP, PEAKS, channel, timePt, ...
    pixelNums, plasmaLams, calLam, xWing, yWing, factor, calPlasmaFit, lineIDs)

% This code loads a plasma shot, fits to two or more lines with known
% wavelengths, and corrects the PIX_SP from the motor calibration.

% addpath('S:\MC_IDS\Matlab Code\Core Codes\lsqcurvefit');
% addpath('S:\MC_IDS\Matlab Code\Calibration\Calibration v2');

% Settings
plotFits = 1;

try
    options = optimsetv61('lsqcurvefit'); % set options for curve fitting
catch
    options = optimset('lsqcurvefit'); % set options for curve fitting
end
f = @singletGauss2D; % function handle

if plotFits
    S = get(0,'ScreenSize');
    fntsz = 18;
    h1 = figure('Visible','on','Name','Individual Fits','Position',...
        [S(3)/6, S(4)/8, 2*S(3)/3, 3*S(4)/4], 'Color', [1 1 1]);
end

% Load Data

% cd('T:\PhantomMovies');
data = importdata(['shot' int2str(shotPlas) '.mat']); % [counts] (time x wavelength space x channel space)
data = squeeze(data(timePt, end:-1:1, end:-1:1));
data = cast(data, 'double');

n_peaks = length(pixelNums);

options.TolFun = 1e-8 * max(data(:)); % set tolerance for curve fitting

centers = NaN*zeros(1, n_peaks);

% Error Analysis:
dCenters = centers;
dp = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]; % fractional increment of 'p' for numerical derivatives

pars = PEAKS(find(PEAKS(:,1) == channel), :);

% Fill in guesses which do not change for each line

guess = [NaN, pars(2), NaN, pars(4), pars(5), NaN];

xBound = round(pars(2)) - xWing : round(pars(2)) + xWing;

for n = 1:n_peaks
    
    % Find approximate 'y' center
    
%     guess(3) = pars(3) - factor * (calLam - plasmaLams(n)) / PIX_SP(find(PEAKS(:,1) == channel));
    guess(3) = pixelNums(n);

    % Make Grid
    
    yBound = round(guess(3)) - yWing : round(guess(3)) + yWing;
    
    [X, Y] = meshgrid(xBound, yBound);
    Z = data(yBound, xBound);
    
    % Refine Initial Guesses
    
    amp = max(Z(:));
    guess(1) = 6 * pars(4) * pars(5) * amp; % volume
    
    guess(6) = min(Z(:)); % offset
    
    % Set Fit Limits
    
    lb = [0.4*guess(1), guess(2)-0.5, guess(3)-1, 0.8*guess(4), 0.9*guess(5), 0];
    ub = [1.8*guess(1), guess(2)+0.5, guess(3)+1, 1.2*guess(4), 3*guess(5), 0.5*amp];
    
    % Reshape data and grid

    x(:, 1) = X(:);
    x(:, 2) = Y(:);
    z = Z(:);

    % Fit Gaussian ----------------------------------------

    parameters = lsqcurvefit(f, guess, x, z, lb, ub, options);
    centers(n) = parameters(3); % save center in wavelength space only
    
    [~,~,dPar,~,~,~,~] = lm(@singletGauss2DLM, parameters, x, z, .001, dp);
    dCenters(n) = dPar(3);
    
    if plotFits
        
        % Plot Data
        
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
        title('Data');
        
        % Reconstruct Fits

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

        ax(2) = axes('Parent', h1, 'Position', [.45 .1 .2 .8], 'FontSize', fntsz);
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
        title('Fit');

        % Link Z Axis

        zmin = min(Z(:));
        zmax = max(Z(:));
        zfmin = min(Zf(:));
        zfmax = max(Zf(:));
        set(ax(1), 'ZLim', [min([zmin zfmin]), max([zmax zfmax])]);
        set(ax(2), 'ZLim', [min([zmin zfmin]), max([zmax zfmax])]);

        pause(1);
        
        if calPlasmaFit.save
            fig_save = getframe(h1);
            [Xfig, mapfig] = frame2im(fig_save);
            imwrite(Xfig, [calPlasmaFit.file num2str(n) '.png']);
        end

    end
    
end

%% Calculate Slope and Adjust PIX_SP

p = polyfit(centers, plasmaLams, 1);

frac = p(1) / PIX_SP(find(PEAKS(:,1) == channel)); % multiplicative factor

PIX_SP = frac * PIX_SP; % adjust PIX_SP
    
%% Plotting
    
if plotFits
    
    h1 = figure('Visible','on','Name','Fit to Plasma Lines','Position',...
        [S(3)/4, S(4)/4, S(3)/2, S(4)/2], 'Color', [1 1 1]);
    
    ax(1) = axes('Parent', h1, 'Position', [.15 .1 .7 .8], 'FontSize', fntsz);
    herrorbar(centers, 1e9 * plasmaLams, dCenters, 'or');
    hold on
    plot(centers, 1e9 * (p(1)*centers + p(2)), '--b');
    grid on;
    for n = 1:length(plasmaLams)
        th = text(centers(n) - 2, 1e9 * (plasmaLams(n) + 0.15 * (plasmaLams(end) - plasmaLams(1))),...
            [{lineIDs{n}; num2str(1e9 * plasmaLams(n))}]);
        set(th, 'FontSize', fntsz -2);
    end
    th = text(centers(2), 1e9 * (plasmaLams(2) - 0.2 * (plasmaLams(end) - plasmaLams(1))),...
        ['PIX SP = ' num2str(p(1)) ' [m/pix]']);
    set(th, 'FontSize', fntsz -2);
    
    title('Measured Dispersion');
    xlabel('Pixels');
    ylabel('Wavelength [nm]');
    
end

end
    
    
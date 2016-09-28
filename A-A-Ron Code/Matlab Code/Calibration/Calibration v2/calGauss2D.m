function[PEAKS, REL_INT, par, fits, guesses] = calGauss2D(peaks, data, brightWing, xWing, indivFit)

% This function loops through each channel and fits a 2D Gaussian, then
% returns the parameters.
%
% par(:,1) = "volume" (the function is a normal distribution)
% par(:,2) = x0
% par(:,3) = y0
% par(:,4) = sigx
% par(:,5) = sigy
% par(:,6) = offset
addpath('T:\RChandra\A-A-Ron Code\Matlab Code\Core Codes\lsqcurvefit');

plotFits = 0;

if plotFits
    S = get(0,'ScreenSize');
    fntsz = 24;
    h1 = figure('Visible','on','Name','Individual Fit','Position',...
        [S(3)/10, S(4)/8, 4*S(3)/5, 3*S(4)/4], 'Color', [1 1 1]);
end

par = NaN*zeros(size(peaks, 1), 6); % allocate fit parameter array
guesses = par; % initial guess array

try
    options = optimsetv61('lsqcurvefit'); % set options for curve fitting
catch
    options = optimset('lsqcurvefit'); % set options for curve fitting
end
options.TolFun = 1e-8 * max(data(:)); % set tolerance for curve fitting
f = @singletGauss2D; % function handle

for n = 1:size(peaks, 1) % loop through channels
    % set grid domain size and create mesh
    
    xBound = round(peaks(n, 2)) - xWing : round(peaks(n, 2)) + xWing;
    yBound = round(peaks(n, 3)) - brightWing : round(peaks(n, 3)) + brightWing;
    
    [X, Y] = meshgrid(xBound, yBound);
    
%     disp(['peaks(' num2str(n) ', 2) = ' num2str(peaks(n, 2))]);
%     disp(['xBound = ' num2str(xBound)]);

    Z = data(yBound, xBound); % Select subset of data for mesh

    if plotFits
        clf; % Plot Raw Data ----------------------------
        ax(1) = axes('Parent', h1, 'Position', [.1 .1 .2 .8], 'FontSize', fntsz);
        h3 = surf(X, Y, Z);
        hold on;
        shading interp;
        colormap jet;
        colorbar;
        view([0 90]);
        set(gca, 'FontSize', fntsz);
        set(gca, 'XLim', [xBound(1) xBound(end)], 'YLim', [yBound(1) yBound(end)]);
        
        xlabel('Real Space');
        ylabel('Wavelength Space');
        title(['Raw Data, Channel ' num2str(peaks(n, 1))]);
    end
    
    % Initial Guesses
    
    sigx = 1.5;
    sigy = 1.0;
    vol = 6 * sigx * sigy * max(Z(:)); % box 'sigma_x' x 'sigma_y' x maximum data point
   
%     [Zmaxy, Zmaxx] = find(Z == max(max(Z))); % find the max point within the grid
%     Zmaxx = round(peaks(n, 2)) - xWing + Zmaxx - 1; % add ZmaxX to the minimum X in the grid
%     Zmaxy = round(peaks(n, 3)) - brightWing + Zmaxy - 1;
    %in theory, this should give the coordinate of the max Z in 
    %terms of the actual cine array.
    
    guess = [vol, peaks(n, 2), peaks(n, 3), sigx, sigy, min(Z(:))];
%     guess = [vol, Zmaxx, Zmaxy , sigx, sigy, min(Z(:))];
    lb = [0.3*vol, guess(2)-0.5, guess(3)-0.5, 0.6*guess(4), 0.6*guess(5), 0]; % lower bound
    ub = [1.8*vol, guess(2)+0.5, guess(3)+0.5, 1.4*guess(4), 1.4*guess(5), 0.5*max(Z(:))]; % upper bound
    
    % Reshape data and grid
    
    x(:, 1) = X(:);
    x(:, 2) = Y(:);
    z = Z(:);

    % Fit Gaussian ----------------------------------------

    [par(n, :), resnorm, residual, exitflag] = ...
        lsqcurvefit(f, guess, x, z, lb, ub, options);
    
    % create fine mesh --------------------------------
    
    nfx = 100;
    nfy = 200;
    xBoundf = linspace(xBound(1), xBound(end), nfx);
    yBoundf = linspace(yBound(1), yBound(end), nfy);
    [Xf, Yf] = meshgrid(xBoundf, yBoundf);

    % Reshape fine mesh for execution by Gaussian function

    xf(:, 1) = Xf(:);
    xf(:, 2) = Yf(:);

    zf = singletGauss2D(par(n, :), xf); % calculate fit on fine mesh

    Zf = reshape(zf, size(Xf, 1), size(Xf, 2)); % reshape into 2D image
    
    if plotFits
        
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
        title(['Gaussian Fit, Channel ' num2str(peaks(n, 1))]);
        
        % Link Z Axis
        
        zmin = min(Z(:));
        zmax = max(Z(:));
        zfmin = min(Zf(:));
        zfmax = max(Zf(:));
        set(ax(1), 'ZLim', [min([zmin zfmin]), max([zmax zfmax])]);
        set(ax(2), 'ZLim', [min([zmin zfmin]), max([zmax zfmax])]);
%         linkprop(ax, 'ZLim');
        
        % Plot Y 2D cross section data and fit ---------------
        
        ax(3) = axes('Parent', h1, 'Position', [.68 .55 .27 .35], 'FontSize', fntsz);
        h5 = plot(yBound, Z(:, ceil(size(Z, 2)/2)), 'or');
        set(h5, 'MarkerFaceColor', 'r');
        hold on;
        h6 = plot(yBoundf, Zf(:, ceil(size(Zf, 2)/2)), '-b');
        set(h6, 'LineWidth', 2);
        set(gca, 'XTickLabel', [], 'YTickLabel', []);
        set(gca, 'XLim', [yBound(1) yBound(end)]);
        grid on;
        title('\lambda Cross Section (y)');
        
        % Plot X 2D cross section data and fit ---------------
        
        ax(4) = axes('Parent', h1, 'Position', [.68 .1 .2 .35], 'FontSize', fntsz);
        h7 = plot(xBound, Z(ceil(size(Z, 1)/2), :), 'or');
        set(h7, 'MarkerFaceColor', 'r');
        hold on;
        h8 = plot(xBoundf, Zf(ceil(size(Zf, 1)/2), :), '-b');
        set(h8, 'LineWidth', 2);
        set(gca, 'XTickLabel', [], 'YTickLabel', []);
        grid on;
        title('Spatial (x)');
        
        % Text -------------
        text(1.08, 0.7, ['\sigma_y = ' num2str(par(n, 5), 4)], 'Units', 'normalized', 'FontSize', fntsz);
        text(1.08, 0.4, ['\sigma_x = ' num2str(par(n, 4), 4)], 'Units', 'normalized', 'FontSize', fntsz);
        
        pause(0.1);
        
        % SAVE
        if and(indivFit.save, indivFit.chan == peaks(n, 1))
            fig_save = getframe(h1);
            [Xfig, mapfig] = frame2im(fig_save);
            imwrite(Xfig, [indivFit.file '.png']);
        end
        
    end % End plotting fits
    
    guesses(n, :) = guess;
    fits.Xf(:, :, n) = Xf; % stack all grids for plotting in main routine
    fits.Yf(:, :, n) = Yf;
    fits.Zf(:, :, n) = Zf;

end % End loop over channels

% Fill in parameters for PEAKS ------------------

PEAKS(:, 1) = peaks(:, 1); % channel numbers
PEAKS(:, 2:5) = par(:, 2:5); % channel centers, x0 and y0, and standard deviations, sigx and sigy

% Convert 'volume' to 'amplitude' to 'area' in y direction only for REL_INT

amp = par(:, 1) ./ (2*pi * par(:, 4) .* par(:, 5)); % volume to amplitude

area = sqrt(2*pi) * par(:, 5) .* amp; % ampitude to area of 1D Gaussian (y only)

REL_INT = mean(area) ./ area; % actually INVERSE relative intensity
    
end
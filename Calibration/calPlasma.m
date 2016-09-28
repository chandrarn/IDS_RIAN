function PIX_SP = calPlasma(shotPlas, PIX_SP, PEAKS, channel, channel2, timePt, ...
    pixelNums, plasmaLams, calLam, xWing, yWing, factor)

% This code loads a plasma shot, fits to two or more lines with known
% wavelengths, and corrects the PIX_SP from the motor calibration.

% Settings
plotFits = 1;

options = optimsetv61('lsqcurvefit'); % set options for curve fitting
f = @singletGauss2D; % function handle

if plotFits
    S = get(0,'ScreenSize');
    fntsz = 18;
    h1 = figure('Visible','on','Name','Individual Fits','Position',...
        [S(3)/6, S(4)/8, 2*S(3)/3, 3*S(4)/4], 'Color', [1 1 1]);
end

% Load Data

cd('T:\PhantomMovies');
data = importdata(['shot' int2str(shotPlas) '.mat']); % [counts] (time x wavelength space x channel space)
timePt
data = squeeze(data(timePt, end:-1:1, end:-1:1));
data = cast(data, 'double');
% figure(9)
% surf(data);
% view([0 90]);
% shading interp;

n_peaks = length(pixelNums);

options.TolFun = 1e-8 * max(data(:)); % set tolerance for curve fitting

centers = NaN*zeros(1, n_peaks);

for k = 1:2 % loop over two channels
    if k == 2
        channel = channel2;
    end
    
    pars = PEAKS(find(PEAKS(:,1) == channel), :);
    
    % Fill in guesses which do not change for each line
    
    guess = [NaN, pars(2), NaN, pars(4), pars(5), NaN];
    
    xBound = round(pars(2)) - xWing : round(pars(2)) + xWing;
    
    for n = 1:n_peaks
        
        % Find approximate 'y' center
        
%         guess(3) = pars(3) - factor * (calLam - plasmaLams(n)) / PIX_SP(find(PEAKS(:,1) == channel));
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
        
        lb = [0.4*guess(1), guess(2)-1.5, guess(3)-1, 0.8*guess(4), 0.9*guess(5), 0];
        ub = [1.8*guess(1), guess(2)+1.5, guess(3)+1, 1.2*guess(4), 3*guess(5), 0.5*amp];
        
        % Reshape data and grid
        
        x(:, 1) = X(:);
        x(:, 2) = Y(:);
        z = Z(:);
        
        % Fit Gaussian ----------------------------------------

        parameters = lsqcurvefit(f, guess, x, z, lb, ub, options);

        centers(n) = parameters(3) % save center in wavelength space only
        
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
            
        end
        
    end
    
    %% Calculate Slope and Adjust PIX_SP
    
    p = polyfit(centers, plasmaLams, 1);
    centers
    PIX_SP(find(PEAKS(:,1) == channel))
    frac = p(1) / PIX_SP(find(PEAKS(:,1) == channel)) % multiplicative factor
    
    % Find index of last fiber in toroidal sub bundle
    lastInd = find(PEAKS(:,1) <= 36, 1, 'last');
    firstInd = find(PEAKS(:,1) >= 37, 1, 'first');
    
    if k == 1
        PIX_SP(1:lastInd) = frac * PIX_SP(1:lastInd); % adjust PIX_SP for toroidal fiber
    else
        PIX_SP(firstInd:end) = frac * PIX_SP(firstInd:end); % adjust PIX_SP for poloidal fiber
    end
end

%% Plotting

if plotFits
    
    h1 = figure('Visible','on','Name','Fit to Plasma Lines','Position',...
        [S(3)/4, S(4)/4, S(3)/2, S(4)/2], 'Color', [1 1 1]);
    
    ax(1) = axes('Parent', h1, 'Position', [.1 .1 .8 .8], 'FontSize', fntsz);
    plot(centers, plasmaLams, 'or');
    hold on
    plot(centers, p(1)*centers + p(2), '--b');
    grid on;
    
    title('Measured Dispersion');
    xlabel('Pixels');
    ylabel('Wavelength [m]');
    
end

end


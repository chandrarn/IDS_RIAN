function[fit_par, dPar, bounds, stddev, param, guesses] = gaussFit_2D(data, param, options)

% Revamped 6-4-13 by ACH to handle 2D fitting

addpath('T:\RChandra\A-A-Ron Code\General Matlab');
addpath('T:\RChandra\A-A-Ron Code\General Matlab\LM Method');

[n_time, n_pix, n_spatial] = size(data);

weight = 0;

fit_par = NaN*zeros(n_time, param.n_chan, 6); % preallocate fit parameters
guesses = fit_par; % preallocate array of initial guesses
bounds = NaN*zeros(n_time, param.n_chan, 4); % preallocate grid bounds array

n_pts = (2*param.xWing + 1) * (2*param.yWing + 1);
blank1 = NaN*ones(1, 6); % overwrite parameters if found to be bad

if param.calcError
    blank2 = NaN*ones(1, n_pts); % overwrite stddev data if bad
    stddev = NaN*zeros(n_time, param.n_chan, n_pts); % preallocate standard deviation
                                                     % array
    dPar = fit_par; % uncertainty in fit parameters

    dp = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]; % fractional increment of 'p' for numerical derivatives
else
    stddev = NaN; dPar = NaN;
end

f1 = @singletGauss2D; % function handle for singlet gaussian function
f2 = @singletGauss2DLM; % LM code must have slightly different form
                    % of function
c = NaN; % not sure if necessary...
% dxData = []; % These could be vectors of known uncertainties in x (pixels)
% dyData = []; %     "                  "               "         y (data values)

% param.Center(52:end) = param.Center(52:end) +1.3;
for n = 1:n_time
    
    for m = 1:param.n_chan%This can be modified to not fit junk fiber data
        
        % Display progress
        if (mod(m,2)==0)
            display(['TIME PT: ' num2str(n) '/' num2str(n_time) ', CHANNEL: ' num2str(m) '/' num2str(param.n_chan)]);
        end
        
        % Make Zero-Velocity Bounds ---------------------------------------
        
        x0 = round(param.peaks(m, 2)-1);
        y0 = round(param.Center(m));% manual offset can go here (I think)
        
        
        xBound = x0 - param.xWing : x0 + param.xWing;
        yBound = y0 - param.yWing : y0 + param.yWing;
        
        % Reset Bounds around Brightest Point------------------------------
%         
%         % Make a tighter bound to avoid resetting the domain too far from
%         % the known zero-velocity line
%         miniBound = yBound(3) : yBound(end - 2);
%         
%         binned = sum(squeeze(data(n, miniBound, xBound)), 2); % bin in real space to give vector in wavelength space
%         
%         [dummy, yMax] = max(binned); % index of maximum value in zero-velocity domain
%         
%         y1 = y0 - param.yWing + yMax + 1; % find y coordinate of brightest point
%         yBound = y1 - param.yWing : y1 + param.yWing; % reset y domain around brightest point
        
        % Save Bounds for Plotting Later ----------------------------------
        
        bounds(n, m, :) = [xBound(1), xBound(end), yBound(1), yBound(end)];
        
        % Make Grid -------------------------------------------------------
        
        [X, Y] = meshgrid(xBound, yBound);
        Z = squeeze(data(n, yBound, xBound));
%         
%         if n >30
%             surf(X,Y,Z); shading interp;
%             view(2);
%             TEMP=inputdlg('TEST: ');
%         end
        
        
        % Initial Guesses -------------------------------------------------
        
        amp = max(Z(:)) - min(Z(:)); % brightest pixel in domain - dimmest

        if amp > param.ampThresh % difference between any two pixels is above threshold

            sigx = param.peaks(m, 4);
            sigy = param.peaks(m, 5);
            vol = 6 * sigx * sigy * amp;

            guess = [vol, param.peaks(m, 2), y0, sigx, sigy, min(Z(:))];
            
            guesses(n, m, :) = guess; % save for post processing

            lb = [0, guess(2)-param.xTol, guess(3)-param.yTol, 0.8*guess(4), 0.9*guess(5), 0];
            ub = [3*vol, guess(2)+param.xTol, guess(3)+param.yTol, 1.2*guess(4), 2.5*guess(5), max(Z(:))];

            % Reshape Data and Grid -------------------------------------------

            x(:, 1) = X(:);
            x(:, 2) = Y(:);
            z = Z(:);

            % Curve Fit -------------------------------------------------------

            [fit_par(n, m, :), resnorm, residual, exitflag] = ...
                lsqcurvefit(f1, guess, x, z, lb, ub, options);

            % Exclude Bad Data ------------------------------------------------

            % Find Amplitude
            amp = fit_par(n, m, 1) / (2*pi * fit_par(n, m, 4) * fit_par(n, m, 5));

            if or(exitflag == 0, amp < param.ampThresh)
                if(and(m>=8,m<=27))
                %disp(['MaxEval reached or Amp<Thresh @: Timepoint: ' num2str(n) ' , Channel: ' num2str(m) ' , Amp: ' num2str(amp) ' , Flag: ' num2str(exitflag) ]);
                end
                fit_par(n, m, :) = blank1;
                if param.calcError

                    dPar(n, m, :) = blank1;
                    stddev(n, m, :) = blank2;
                end
            elseif param.calcError
%                 % Find error estimates using modified LM code
%                 [X2, dPar(n, m, :), stddev(n, m, :), corr, R_sq] = ...
%                     lm_sigma(f2, fit_par(n, m, :), x, z, dp, c, weight);
                   fit_par(n,m,1)=fit_par(n,m,1)-3;
                  [p_fit,Chi_sq,dPar(n,m,:),stddev(n,m,:),corr,R2,cvg_hst] = ...
                    lm(@singletGauss2DLM, fit_par(n,m,:), x, z, .001,dp);%, p_min,p_max,0)
                    
            end

        else % Do not attempt to fit
            if(and(m>=8,m<=27))
            %disp(['Amp <= Threshold @: Timepoint: ' num2str(n) ' , Channel: ' num2str(m) ]);
            end
            fit_par(n, m, :) = blank1;
            if param.calcError
                dPar(n, m, :) = blank1;
                stddev(n, m, :) = blank2;
            end
        end % Attempt the fit at all
        
    end % Channel Loop
    %assignin('base','fit_par',fit_par);
end % Time Loop
end
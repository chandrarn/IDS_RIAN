function[fit_par, dPar, bounds, stddev, param, guesses] = gaussFit_2D(data, param, options, nn)
global homePath
% Revamped 6-4-13 by ACH to handle 2D fitting
% fit_par = [vol, x0, y0, sigx, sigy, min(Z(:))];
%addpath('T:\RChandra\A-A-Ron Code\General Matlab');
%addpath([homePath '\A-A-Ron Code\Matlab Code\Core Fitting Codes\lsqcurvefit\']);%Git Compliant
%addpath('T:\RChandra\NewCodes\Matlab\lsqcurvefit\');
%addpath('T:\RChandra\NewCodes\Matlab\lM Method\');
% addpath('T:\RChandra\NewCodes\Matlab\');


[n_time, ~, ~] = size(data); % time, wavelength, spatial

% weight = 0;

fit_par = NaN * zeros(n_time, param.n_chan, 6); % preallocate fit parameters
guesses = fit_par; % preallocate array of initial guesses
bounds = NaN * zeros(n_time, param.n_chan, 4); % preallocate grid bounds array

n_pts = (2 * param.xWing + 1) * (2 * param.yWing + 1);
blank1 = NaN * ones(1, 6); % overwrite parameters if found to be bad

if param.calcError
    blank2 = NaN * ones(1, n_pts); % overwrite stddev data if bad
    stddev = NaN * zeros(n_time, param.n_chan, n_pts); % preallocate standard deviation
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
% figure; surf(squeeze(sum(data(:,mean(param.Center(:, nn))-param.yWing:mean(param.Center(:, nn))+param.yWing,:),1))); shading interp;
% view([ 0 90]);set(gca,'ylim',[0,size(data,2)]);
fig =0;
for n = 1:n_time
    
    for m = 1:param.n_chan
        %progress
        if(mod(m,10)==0)
            display(['CHANNEL: ' num2str(m) ' OF: ' num2str(param.n_chan) '. TIME: ' num2str(n) ' OF: ' num2str(n_time)]);
        end
        % Make Zero-Velocity Bounds ---------------------------------------
        
        x0 = round(param.peaks(m, 2));
        y0 = round(param.Center(m, nn));
        
%         % FOR 151217026 2/24/16 ELIMINATE
%         if ~isempty(find(m == [16,22,26]))
%             y0=y0 -.5;
%         end
%         if m == 26
%             y0=y0-.5;
%         end
%         
%         y0=round(y0);
        
        xBound = x0 - param.xWing : x0 + param.xWing;
        if nn == 2
            yBound = y0 - param.yWing-4 : y0 + param.yWing+1; %ELIMINATE -3 23/2/16
        elseif nn == 3
            yBound = y0 - param.yWing-4 : y0 + param.yWing-1;
        elseif nn == 4
            yBound = y0 - param.yWing +0: y0 + param.yWing-4;
        else
            yBound = y0 - param.yWing -1: y0 + param.yWing+1;
        end
        yBound = yBound -1*(m>29); % second set of fibers is slightly lower.
        if param.calcError
            n_pts = (2 * param.xWing + 1) * (length(yBound));
            blank1 = NaN * ones(1, 6); % overwrite parameters if found to be bad
            blank2 = NaN * ones(1, n_pts); % overwrite stddev data if bad
            stddev = NaN * zeros(n_time, param.n_chan, n_pts);
        end
        
        if fig ==0 %&& m>29
            [X,Y] = meshgrid(1:size(data,3),yBound);
            yBound
            figure; surf(X,Y,squeeze(sum(data(:,yBound,:),1))); shading interp;
            view([ 0 90]);set(gca,'ylim',[0,size(data,2)]);
            fig =1;
        end

        
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
        
            
        
        
        
        % Initial Guesses -------------------------------------------------
        
        amp = max(Z(:)) - min(Z(:)); % brightest pixel in domain - dimmest

        if amp > param.ampThresh % difference between any two pixels is above threshold
            
            % test autofinding y0,x0
            [x1,y1] = maXY(Z);
%             [~,I] = max(max(Z,[],2));
            y1 = y1-1 + yBound(1);
%             [~,I] = max(max(Z,[],1));
            x1 = x1-1 + xBound(1);
            if(m==26)
                display(['Error, (x,y) = (' num2str(x1-x0) ',' num2str(y1-y0) ')']);
                amp = amp *1.1; % ELIM 2/24/16, for 151217026
            end
            
            x0=x1;
            y0=y1;
            
            sigx = param.peaks(m, 4);
            sigy = param.peaks(m, 5);
            vol = 2*pi * sigx * sigy * amp;
            % changed 2/24/16, testing using maXY x0,y0 values
            %guess = [vol, x0, y0, sigx, sigy, min(Z(:))];
            guess = [vol, param.peaks(m, 2), y0, sigx, sigy, min(Z(:))];
            % ELIMINATE 21/1/16
            if find(m ==[14:18,22-3,24-3,26-3,57-7,58-7,49-7,50-7])
                guess(1) = guess(1)*1.4;
                guess(6) = guess(6)*2.3;
                guess(5) = guess(5).*1.2;
            end
            % Give the temperature guess a little boost to start
            guess(5) = 1.5 * guess(5);
            
            guesses(n, m, :) = guess; % save for post processing

            lb = [0, guess(2)-param.xTol, guess(3)-param.yTol, 0.8*guess(4), 0.5*guess(5), 0];
            ub = [3*vol, guess(2)+param.xTol, guess(3)+param.yTol, 1.2*guess(4), 4*guess(5), max(Z(:))];
            
            
%             assignin('base','guess',guess);
%             figure;
%             surf(X,Y,Z); shading interp;
%             view([0,90]);
%             TEMP=inputdlg('TEST: ');
%             if TEMP{1} == '1'
%                 myVarList=who;
%                 for indVar = 1:length(myVarList)
%                     assignin('base',myVarList{indVar},eval(myVarList{indVar}));
%                 end
%                 return;
%             end

            % Reshape Data and Grid -------------------------------------------
            x(:, 1) = X(:);
            x(:, 2) = Y(:);
            z = Z(:);

            % Curve Fit -------------------------------------------------------
             try
            [fit_par(n, m, :), resnorm, residual, exitflag] = ...
                lsqcurvefit(f1, double(guess), double(x), double(z), double(lb), double(ub), options);
             catch
                exitflag=0;
                 display(['Lsqcurvefit error: ' num2str(n) ',' num2str(m)]);
             end
            % Exclude Bad Data ------------------------------------------------

            % Find Amplitude
            amp = fit_par(n, m, 1) / (2*pi * fit_par(n, m, 4) * fit_par(n, m, 5));

            if or(exitflag == 0, amp < param.ampThresh)
                if(and(m >= 8,m <= 27))
                disp(['MaxEval reached or Amp<Thresh @: Timepoint: ' num2str(n) ' , Channel: ' num2str(m) ' , Amp: ' num2str(amp) ' , Flag: ' num2str(exitflag) ]);
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
                   fit_par(n, m, 1) = fit_par(n, m, 1) - 3;
                  [p_fit, Chi_sq, dPar(n, m, :), stddev(n, m, :), corr, R2, cvg_hst] = ...
                    lm(@singletGauss2DLM, fit_par(n, m, :), x, z, 0.001, dp);%, p_min,p_max,0)
                    
            end

        else % Do not attempt to fit
%             if(and(m >= 8,m <= 27))
%             %disp(['Amp <= Threshold @: Timepoint: ' num2str(n) ' , Channel: ' num2str(m) ]);
%             end
            fit_par(n, m, :) = blank1;
            if param.calcError
                dPar(n, m, :) = blank1;
                stddev(n, m, :) = blank2;
            end
        end % Attempt the fit at all
        
    end % Channel Loop
end % Time Loop
end
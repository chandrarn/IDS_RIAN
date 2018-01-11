function [PIX_SP,parameters,data] = fitMotor(shot3, shot4, PEAKS, motorSpeed, brightWing, xWing, trimData)
%addpath('S:\MC_IDS\Matlab Code\Core Codes\lsqcurvefit');
% SPECIAL NOTE:
% +1 ADDED TO PEAKS, TIME BOUNDS ADDED TO MOVIE
% Updated to return all fit parameters for all timepoints
plotFits = 0;

% Load Data
%PEAKS(2,:) = PEAKS(2,:)+1;
addpath('T:\PhantomMovies');

% try
     data = importdata(['Shot ' int2str(shot3) '.mat']);
     time3 = data.TimeVector;
     data = data.CineArray;
     % This is necessary for shitty, shitty back compatibility
     data = shiftdim(data,2);
     
     if isempty(trimData)
         trimData1 = 1:size(data,1);
     else
         trimData1=trimData;
     end
%      Unfortunately, the following all must be within the TRY because the
%      oldformat is moronic, and put time as the first column.
% The sizes here sometimes need to be manually changed, IE: avi converter
% may do the indexing differently. note 'shiftdim'
    
     [ n_time,n_pix, n_chan,] = size(data);
     data3 = zeros(length(trimData1), n_pix, n_chan);
     for i = 1:length(trimData1)
         data3(i, 1:n_pix, 1:n_chan) = squeeze(data(trimData1(1) + i - 1,:,:));
     end
%      data3 = shiftdim(data,3);
     data3 = data3(:, end:-1:1, :);%(trimData,:,:);
     time3=time3(trimData1);
% catch
%     
%     data3 = importdata(['shot' int2str(shot3) '.mat']); % [counts] (time x wavelength space x channel space)
%     time3 = importdata(['t' int2str(shot3) '.mat']); % [ms]
%     
%     if isempty(trimData)
%         trimData = 1:size(data3, 1);
%     end
%     data3 = data3(trimData, end:-1:1, end:-1:1);
%     time3 = time3(trimData);
%     clear h1;
% end
[n_time, n_pix, n_chan] = size(data3);
assignin('base','data3',data3);
% THE EXPLICIT TIME TRIMMING IS BECAUSE YOU DONT WANT TO FIT TO THE SECOND
% LINE

if shot4 ~= 0 % using both fibers
%     try
         data = importdata(['Shot ' int2str(shot4) '.mat']);
         time4 = data.TimeVector;
         data = data.CineArray;
         data = shiftdim(data,2);
         if isempty(trimData)
             trimData2 = 1:size(data,1);
         else
             trimData2=trimData;
         end
%          Unfortunately, the following all must be within the TRY because the
%          oldformat is moronic, and put time as the first column.
         [n_time,n_pix, n_chan] = size(data);
         data4 = zeros(length(trimData2), n_pix, n_chan);
         for i = 1:length(trimData2)
             data4(i, 1:n_pix, 1:n_chan) = squeeze(data( trimData2(1) + i - 1,:,:));
         end
         time4=time4(trimData2);
%          data4 = shiftdim(data,3);
         data4 = data4(:, end:-1:1, :);%(trimData,:,:);
%     catch
% 
%         data4 = importdata(['shot' int2str(shot4) '.mat']); % [counts] (time x wavelength space x channel space)
%         time4 = importdata(['t' int2str(shot4) '.mat']); % [ms]
% 
%         if isempty(trimData)
%             trimData = 1:size(data4, 1);
%         end
%         data4 = data4(trimData, end:-1:1, end:-1:1);
%         time4 = time4(trimData);
%         clear h1;
%     end
    % Add movies together

    if length(time3) > length(time4)
                size(time3)
%         size(time4)
%         size(data3)
%         size(data4)
        data = data4 + data3(1:length(time4), :, :);
        time = time4;
        clear data3 data4 time3 time4;
    else
%         size(time3)
%         size(time4)
%         size(data3)
%         size(data4)
        data = data3 + data4(1:length(time3), :, :);
        time = time3;
        clear data3 data4 time3 time4;
    end
else % only using one fiber
    data = data3;
    time = time3(trimData);
    clear data3 time3;
end

%cd('S:\MC_IDS\Matlab Code\Calibration\Calibration v2');

data = cast(data, 'double');
[n_time, n_pix, n_chan] = size(data)

%% Loop over Channels
n_chan = size(PEAKS,1);

options = optimsetv61('lsqcurvefit'); % set options for curve fitting
options.TolFun = 1e-8 * max(data(:)); % set tolerance for curve fitting
f = @singletGauss2D; % function handle

centers = NaN * zeros(n_chan, length(time)); % initialize array to store centers for all time
%h2 = figure;
if plotFits
    S = get(0,'ScreenSize');
    fntsz = 20;
    h1 = figure('Visible','on','Name','Individual Fit','Position',...
        [S(3)/10, S(4)/8, 4*S(3)/5, 3*S(4)/4], 'Color', [1 1 1]);
end

finalfit = 0;

% lm error fitting
%addpath('T:\RChandra\A-A-Ron Code\General Matlab');
%addpath('T:\RChandra\A-A-Ron Code\General Matlab\LM Method');
fit_par = NaN*zeros(length(time),n_chan,2);
stddev = NaN*zeros(length(time),n_chan); % preallocate standard deviation
                                                 % array
dPar = fit_par; % uncertainty in fit parameters

dp = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001];
f1 = @linearFit;
f2 = @linearLM;
c= NaN;

parameters = zeros(n_time,n_chan,6); % Storing parameters

for n = 1:n_chan

    % find real space bound for fitting domain

    xBound = round(PEAKS(n, 2)) - xWing : round(PEAKS(n, 2)) + xWing;

    for m = 1:n_time

        % display progress
         if (mod(m, 10) == 0) 
            display(['CHANNEL: ' num2str(n) '/' num2str(n_chan) ', TIME PT: ' num2str(m) '/' num2str(n_time)]);
         end
        % find wavelength coordinate of starting peak
        [amp, col] = max(max(data(m, :, xBound))); % which column has the max
        [amp, y] = max(squeeze(data(m, :, xBound(col)))); % which row in that column
        

        % find wavelength direction bound for fitting
        %pause(1.5);
        yBound = y - brightWing : y + brightWing;
        while yBound(1) <1 % make sure we arent looking off CCD
            yBound = yBound +1;
        end
        while yBound(end) >n_pix 
            yBound = yBound -1;
        end
        %col = PEAKS(n,2)
         % test quadratic fitting for extra fittness.
         PEAKS(n,5) * amp * 3;
         y;
         PEAKS(n, 5);
         assignin('base', 'nYbounds', yBound);
         assignin('base', 'm', m);
         min(data(m, yBound, round(PEAKS(n,2))));
        guess =[PEAKS(n,5)*amp*3, y, PEAKS(n,5), min(data(m, yBound, round(PEAKS(n, 2))))];
        assignin('base','guess1',guess);
        
        PreParam = lsqcurvefit(@singletFun,guess,yBound,data(m,yBound,round(PEAKS(n,2))));
        
        %origy= y
        y = PreParam(2);
        
        %close all;
        %ax(3) = axes('Parent',h2);
        % create fine mesh --------------------------------

%             [Xf, Yf] = meshgrid(xBound, yBound);
% 
%             % Reshape fine mesh for execution by Gaussian function
% 
%             xf(:, 1) = Xf(:);
%             xf(:, 2) = Yf(:);
        %surf(Xf,Yf,squeeze(data(m,yBound,xBound))); shading interp; view([ 0 90]);
        %figure;plot(yBound,singletFun(PreParam,yBound)); hold on;
        %plot(yBound,data(m,yBound,round(PEAKS(n,2))),'r*');
        
%         hold off;
        
%         assignin('base','PreParam',PreParam);
%         assignin('base','guess',guess);

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
        
        if plotFits %&& mod(m,4) ==0
            clf; % Plot Raw Data ----------------------------
            ax(1) = axes('Parent', h1, 'Position', [.1 .1 .2 .8], 'FontSize', fntsz);
            h3 = surf(X, Y, Z);
            hold on;
            shading interp;
            colormap jet;
            colorbar;
            caxis([ 0 4500]);
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
        
        guess = [vol, PEAKS(n, 2)-.5, y, sigx, sigy, min(Z(:))];
        assignin('base','guess2',guess);

        lb = [0.4*vol, guess(2)-0.5, guess(3)-0.5, 0.8*guess(4), 0.8*guess(5), 0]; % lower bound
        ub = [1.8*vol, guess(2)+0.5, guess(3)+0.5, 1.2*guess(4), 1.2*guess(5), 0.5*amp]; % upper bound
        if lb == ub
            display('FITTING BOUNDS EQUAL, TERMINATING');
        end
        assignin('base','lb',lb);
        assignin('base','ub',ub);
        % Reshape data and grid

        x(:, 1) = X(:);
        x(:, 2) = Y(:);
        z = Z(:);

        % Fit Gaussian ----------------------------------------
        
        parameters(m,n,:) = lsqcurvefit(f, guess, x, z, lb, ub, options);
       % parameters1 = lsqcurvefit(f, [vol, PEAKS(n, 2)-.5, origy, sigx, sigy, min(Z(:))], x, z, lb, ub, options);
%        origfit = parameters1(3);
        if parameters(m,n,3) == finalfit
            display(['NOCHANGECONDITIONREACHED: T ' num2str(m) ', P ' num2str(n)]);
        end
        finalfit = parameters(m,n,3);
        centers(n, m) = parameters(m,n,3); % save center in wavelength space only
        
        if plotFits %&& mod(m,4) ==0
            % create fine mesh --------------------------------
            nfx = 100;
            nfy = 200;
            xBoundf = linspace(xBound(1), xBound(end), nfx);
            yBoundf = linspace(yBound(1), yBound(end), nfy);
            [Xf, Yf] = meshgrid(xBoundf, yBoundf);

            % Reshape fine mesh for execution by Gaussian function

            xf(:, 1) = Xf(:);
            xf(:, 2) = Yf(:);

            zf = singletGauss2D(parameters(m,n,:), xf); % calculate fit on fine mesh

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
    
end % channel loop

assignin('base','centers',centers);

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
Chi_sq = zeros(n_chan,1);
dPar = zeros(n_chan,2);
stddev = zeros(n_chan,length(time));
corr = zeros(n_chan,2,2);
R = zeros(n_chan,1);
for n = 1:n_chan
    % slope estimate nm/s * s/pix = [nm/pix]
%     A(1) = dnm/(centers(n, end) - centers(n, 1));
%     A(2) = nm(1); % offset estimate
%     [estimate] = optimize(centers(n, :)', nm, A, max_evals);
n
    centers(n, :)';
    nm;
    size(centers)
    size(nm)
    p = polyfit(centers(n, :), nm', 1);
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

     %[p_fit(n,:),Chi_sq(n),dPar(n,:),stddev(n,:),corr(n,:,:),R(n)] = ...
      %              lm(f2, guess, x, z, .001,dp);%, p_min,p_max,0)
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








%% Correct PIX_SP by finding the motor speed and the pix/sec
% uses a doublet to find the motor speed vL (v lambda ), then looks at
% the actual data location, finds the line speed in pixels per second,
% calculates PIX_SP in lam/sec, corrects the given values for PIX_SP

%NOTE: previously, was not fitting whole line, just the line connecting the
%start and endpoints. This is not a terrible model, but still, isnt as
%robust as it could be.
function PIX_SP = newCalPlasma(motorCalShot, PIX_SP, binChanMotor, lamMotor,  PEAKS, channel, xWing)
        
        % if the doublet motor lines are from a known shot, add it here for
        % speed of use.
%         if motorCalShot == 2222222
%             vL = .9883;
%             realPIX = 1.1601e-11;
%         else
            addpath('T:\PhantomMovies');
            data = importdata(['shot' int2str(motorCalShot) '.mat']); % [counts] (time x wavelength space x channel space)
            time = importdata(['t' int2str(motorCalShot) '.mat']); % [ms]
            data = cast(data, 'double');
            data = data(:, end:-1:1, end:-1:1);
            [n_time, n_wavelength, n_spatial] = size(data);
            for n = 1:n_time
                 dataB(n, :) = sum(squeeze(data(n, :, binChanMotor-xWing:binChanMotor+xWing)), 2);
            end
            
            %% Fit Gaussians from Motor Shot
            options = optimsetv61('lsqcurvefit'); % set options for curve fitting
            options.TolFun = 1e-8 * max(dataB(:)); % set tolerance for curve fitting
            f = @singletFun; % function handle
            
            
            motorCenters = zeros(n_time, size(binChanMotor, 1));
            xdata = 1:n_wavelength; % 1:96
            for n = 1:n_time % loop over time to find the speed of the line
                ydata = squeeze(dataB(n, :)); % curve for one channel, one time

                [amp, x0] = max(ydata);
                fwhm = length(find(ydata > 0.5 * amp));
                sigma = fwhm / 2.35;
                area = 0.8 * sqrt(2*pi) * sigma * amp; % 0.8 factor arbitrary for best initial guess
                offset = min(ydata);

                guess = [area, x0, sigma, offset];
                lb = [0.3*area, guess(2)-1, 0.2*guess(3), 0]; % lower bound
                ub = [3*area, guess(2)+1, 3*guess(3), 0.5*amp]; % upper bound
                parameters = lsqcurvefit(f, guess, xdata, ydata, lb, ub, options);

                motorCenters(n) = parameters(2);
                if mod(n,100) ==0
                    display(['Time Point: ' num2str(n) '/' num2str(n_time)]);
                end
            end
            
% %             % Plot position vs frame to find the location of the two lines
%             figure('units','normalized','outerposition',[0 0 1 1]); hold on;
             figure; plot(motorCenters);hold on;
             assignin('base','motorCenters',motorCenters);
%             timePoints(1,1) = input('First Frame of first line: ');
%             timePoints(1,2) = input('Second Frame of first line: ');
%             timePoints(2,1) = input('First Frame of second line: ');
%             timePoints(2,2) = input('Second Frame of second line: ');
% % Fit the two lines
%             [P1,S1] = polyfit(time(timePoints(1,1:2))*1e-3,motorCenters(timePoints(1,1:2)),1); 
%             [P2,S2] = polyfit(time(timePoints(2,1:2))*1e-3,motorCenters(timePoints(2,1:2)),1);
%             assignin('base','P1',P1);
%             assignin('base','P2',P2);
%             x = input('test');
% auto detect two lines: move through x axis till you find a
            % run of 30 frames that have a reasonable slope, and then find
            % when it ends, then do that again for the second line, fully
            % fit the two lines, get the two slopes.
            frame0 = 0; % start of first line
            frame1 = 0; % end of first line
            frame2 = 0;
            frame3 = 0;
            %slope range only valid for max speed motor
            %slopeRange = [90,80]; % approximate max and min values that the slope of the line can have
            %slopeRange = [20,10]; % for 12nm/min
            slopeRange = [5,3]; %
            saftey = 25; % saftey factor, the first few frames of the line are generaly poor.
            %Note: can change '50' if the lines are really short;
            for i = 1:size(motorCenters,1)-50
                [P,S] = polyfit(time(1:50)*1e-3,motorCenters(0+i:49+i)',1);
                Q(i)=P(1);
                if (P(1) < slopeRange(1) && P(1) > slopeRange(2)) % if we see the line
                    if frame0 == 0 % first time we see it
                        frame0 = i+saftey% move it forward a little
                    elseif frame1 ~= 0 && frame2 == 0 % second time we're seeing it
                        frame2 = i+saftey
                    end
                elseif ~(P(1) < slopeRange(1) && P(1) > slopeRange(2)) % if we loose the line
                    if frame0 ~= 0 && frame1 == 0 && frame2 == 0% we've seen it once
                        frame1 = i+49-saftey% then we just lost it
                        [P1,S1] = polyfit(time(frame0:frame1)*1e-3,motorCenters(frame0:frame1)',1);
                    elseif frame2 ~= 0 && frame3 == 0 % we've seen it twice
                        frame3 = i+49-saftey; % then we've just lost it again
                        [P2,S2] = polyfit(time(frame2:frame3)*1e-3,motorCenters(frame2:frame3)',1);
                    end
                end
            end
                    
            assignin('base','Q',Q);
            assignin('base','P',[P1;P2]);
            plot(frame0,motorCenters(frame0),'r*',frame1,motorCenters(frame1),'r*',frame2,motorCenters(frame2),'r*',frame3,motorCenters(frame3),'r*');
            
            % loop over y axis to find average distance x between lines
            % line will be about 100 frames long ( assuming same framerate
            % is used as shot 2222222 ) and we dont like the data at the
            % very bottom of the screen ( hense +50 ). Solve the two lines
            % to find the x points, from y = 50 to y = 150.
            for i = 0:100
                dt(i+1) = ((i+50 -P2(2))/P2(1)) - ((i+50 -P1(2))/P1(1));
            end
            dt = mean(dt);
            dL = (lamMotor(2)-lamMotor(1)); % change in Lambda in nm
            vL = dL/dt % lambda per second dL/dt = MOTOR SPEED
            vP = mean([P1(1),P2(1)]);
            
%         end
        
        
          % Old version: calculated pixels per second with a seperate movie
          % Although the other movie was at a lambda closer to the data
          % collection tuning, we've found that PIX SP doesnt change that
          % terribly much at different spectrometer settings.
          % 
%         if plasmaCalShot == 3333333
%             vP = 85.8998;
%         else
%             addpath('T:\PhantomMovies');
%             data = importdata(['shot' int2str(plasmaCalShot) '.mat']); % [counts] (time x wavelength space x channel space)
%             time = importdata(['t' int2str(plasmaCalShot) '.mat']); % [ms]
%             data = cast(data, 'double');
%             data = data(:, end:-1:1, end:-1:1);
%             [n_time, n_wavelength, n_spatial] = size(data);
%             for n = 1:n_time
%                  dataB(n, :) = sum(squeeze(data(n, :, binChanPlas-xWing:binChanPlas+xWing)), 2);
%             end
%             
%             %% Fit Gaussians from Motor Shot
%             options = optimsetv61('lsqcurvefit'); % set options for curve fitting
%             options.TolFun = 1e-8 * max(dataB(:)); % set tolerance for curve fitting
%             f = @singletFun; % function handle
%             
%             
%             motorCenters = zeros(n_time, size(binChanPlas, 1));
%             xdata = 1:n_wavelength; % 1:96
%             for n = 1:n_time % loop over time to find the speed of the line
%                 ydata = squeeze(dataB(n, :)); % curve for one channel, one time
% 
%                 [amp, x0] = max(ydata);
%                 fwhm = length(find(ydata > 0.5 * amp));
%                 sigma = fwhm / 2.35;
%                 area = 0.8 * sqrt(2*pi) * sigma * amp; % 0.8 factor arbitrary for best initial guess
%                 offset = min(ydata);
% 
%                 guess = [area, x0, sigma, offset];
%                 lb = [0.3*area, guess(2)-1, 0.2*guess(3), 0]; % lower bound
%                 ub = [3*area, guess(2)+1, 3*guess(3), 0.5*amp]; % upper bound
%                 parameters = lsqcurvefit(f, guess, xdata, ydata, lb, ub, options);
% 
%                 motorCenters(n) = parameters(2);
%             end
%             
%             % Plot position vs frame to find the location of the two lines
%             figure('units','normalized','outerposition',[0 0 1 1]); hold on;
%             plot(motorCenters);
%             timePoints(1,1) = input('First Frame of first line: ');
%             timePoints(1,2) = input('Second Frame of first line: ');
% 
%             
%             % Fit the two lines
%             [P1,S1] = polyfit(time(timePoints(1,1):timePoints(1,2))*1e-3,motorCenters(timePoints(1,1):timePoints(1,2)),1); 
%             
%             vP = P1(1) % pixels per second
%         end
        
        
         % using pixel speed and motor speed, can computer PIX SP
        dLdP = vL/vP;

        NewPIX_SP = dLdP * 1e-9
        

        frac = NewPIX_SP / PIX_SP(find(PEAKS(:,1) == channel)) % multiplicative factor
        PIX_SP = frac*PIX_SP; % return scaled PIX_SP
end
       
                
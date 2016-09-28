%% Remove Slant-step offset(signal, t, shot, Signal Name)
function [signal, t] = Subtract_Step_Slat_Offset(signal, t, shot, Sig_Name)
% Removes a slant-step from a signal by presuming a steady-state signal before and
% after the injectors are powered, and a linear increase or decrease in the
% signal during the period the injectors are active.  
% All signals on HIT-SI that are integrated (DAFI or Numeric) see a
% step-slant offset of some size.  Likely a rise in the ground during the
% shot or something similar.
% So first I'm  going to integrate the signal to do the calculations, but
% then I'm going to do the subtraction on the original un-integrated signal
% that was passed in.  This means that I'm going to subtract off a top-hat
% funtion instead of a ramp function.
% written by JSW
    
global OUTPUT

    [a status] = mdsopen('hitsi',shot);                                    % Open the shot
        active_x   = mdsvalue('\PSI_INJ_X');                               % Load a signal that defines the start and stop of the shot (flux,I_pri, a voltage divider on the back on an SPA, or SPA_ENABLE?
        active_x_t = mdsvalue('dim_of(\PSI_INJ_X)');                       % PSI_INJ is a good signal because it turns on and off FAST and is only on during the shot
        active_y   = mdsvalue('\PSI_INJ_Y');                               % Load the Y side as well, PSI_INJ_X is in units of Wb
        active_y_t = mdsvalue('dim_of(\PSI_INJ_Y)');                       % And the time basis
        if (max(abs(active_x)) < .1e-3 ) && (max(abs(active_y)) < .1e-3)   % If the shot is a voltage coil only shot (flux less than .1mWb) then use a substitute signal: Voltage demand
            active_x     = mdsvalue('\X_VOLT_DEM');                        % Load a different signal that defines the start of the shot (flux,I_pri, a voltage divider on the back on an SPA, ro SPA_ENABLE?
            active_x_t   = mdsvalue('dim_of(\X_VOLT_DEM)');                % PSI_INJ is a good signal because it turns on and off FAST and is only on during the shot
            active_y     = mdsvalue('\Y_VOLT_DEM');                        % Load the Y side as well
            active_y_t   = mdsvalue('dim_of(\Y_VOLT_DEM)');                % And the time basis
        end                                                                %         
    mdsclose;                                                              % Close the shot

    % Window the shot in time:
    start_ind_x = find(abs(active_x)>5e-5, 1, 'first');                    % Returns an index for the start of the shot: first time PSI_INJ_X is above .05 mWB
    start_ind_y = find(abs(active_y)>5e-5, 1, 'first');                    % Returns an index for the start of the shot: first time PSI_INJ_Y is above .05 mWB
    stop_ind_x  = find(abs(active_x)>5e-5, 1, 'last');                     % Returns an index for the stop of the shot: last time PSI_INJ_X is above .05 mWb
    stop_ind_y  = find(abs(active_y)>5e-5, 1, 'last');                     % Returns an index for the stop of the shot: last time PSI_INJ_Y is above .05 mWb  
    start_t_x   = active_x_t(start_ind_x);                                 % Inline equations slow matlab down, that's why this is broken apart
    start_t_y   = active_y_t(start_ind_y);                                 % Convert the start indicies into time values
    stop_t_x    = active_x_t(stop_ind_x);                                  % Convert the stop indicies into time values
    stop_t_y    = active_y_t(stop_ind_y);                                  % 
    start_t     = min([start_t_x, start_t_y]);                             % Find the earliest start time to use as the shot start time
    stop_t      = max([stop_t_x, stop_t_y]);                               % Find the latest stop time to use as the shot end time

    if isempty(start_t) || isempty(stop_t)                                 % Check to see if there was a plasma ever or if it was a vcoil only shot
        start_t = 0;                                                       % If it's a vcoil only shot then we need to pick a shot length
        stop_t  = 10e-3;                                                   %  instead of trying to determine the shot length (could use SPA_ENABLE)
    end                                                                    % 
    
    % Calculate the step-slant offset
    [junk,start_ind] = min(abs(t - start_t));                              % Find the index for the probe signal that corresponds to the shot starting time
    [junk,stop1_ind] = min(abs(t - stop_t));                               % Find the index for the probe signal that corresponds to the shot ending time
    [junk,stop2_ind] = min(abs(t - (stop_t + 2e-3)));                      % A second stop time that gives the shot 2ms to settle back down to zero
    [sig_int, t_int] = Integrate(signal,t);                                % Numerically integrate the raw-but-filtered probe signal        
    y1 = polyfit(t(1:start_ind), sig_int(1:start_ind),0);                  % Fit a horizontal line to the first part (t<0) 
    y2 = polyfit(t(stop2_ind:end), sig_int(stop2_ind:end),0);              % Fit a horizontal line to the last part (t>end_of_shot+2ms)
    
    % Create the step-slant subtraction signal
    dh = y2-y1;                                                            % This is how much of a step occurs though the shot
    dt = stop_t-start_t;                                                   % Find the run, as in rise/run, to determine the slope
    h  = dh/dt;                                                            % Integrating h*dt is the area or dh that the step-slant would make it up to
    y_sub = [zeros(1,start_ind-1), h*ones(1,stop1_ind-start_ind+1),...     %
               zeros(1,length(signal)-stop1_ind)]';                        %
    
    if OUTPUT > 1
        got_it = 0;                                                        % Whether or not the code has found an unused Figuer(20_)figure
        fig_num = 0;                                                       % Variable to incrementing Figure(200)
        while ~got_it                                                      % Keep running until an empty figure is found
            try                                                            % try/catch because testing for an empty figure requires get() to bomb out
                junk=get(200+fig_num);                                     % This call will bomb out if the figure doesn't exist
                fig_num = fig_num + 1;                                     % 
            catch                                                          %
                got_it = 1;                                                % If get( ) bombs out then that's the open figure I want to write to
            end                                                            % 
        end                                                                % 
        [sss_int, t_sss]   = Integrate(signal - y_sub, t);                 % 
        [y_sub_int, t_sss] = Integrate(y_sub, t);                          % 
        figure(200+fig_num)                                                % Opening the next new figure
            ax(1) = subplot(2,1,1);                                        % 
                plot(t,signal), hold all, grid on                          % 
                plot(t,(signal - y_sub)), xlabel('Time (s)'), ylabel('Signal (V?)')
                plot(t, y_sub)                                             % Plot the step offset
                title([Sig_Name ': Filtered Signal Before and After Step Slant Subtraction, Shot: ' num2str(shot)]);
                legend('Filtered','Step-Slant Subtracted','Step');         % 
            ax(2) = subplot(2,1,2);                                        % 
                plot(t,sig_int), hold all, grid on                         %
                plot(t,sss_int), xlabel('Time (s)'), ylabel('Int(Signal) (Vs?)');
                plot(t, y_sub_int)                                         % Plot the slant offset
                title([Sig_Name ': Integrated Signal, Before and After Step Slant Subtraction, Shot: ' num2str(shot)]);
                legend('Filtered','Step-Slant Subtracted','Step-Slant');   % 
            linkaxes(ax,'x')
    end
     
    % RETURN: [signal, t]
    if 1                                                                   % If the step-slant is substantial then
        signal = signal - y_sub;                                           %   subtract the offset
    end                                                                    % Otherwise the step-slat is small and just return the original signal

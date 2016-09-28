%% function [I, t, tag, ProbeProfile, AmperProfile, ZerosProfile] = Calc_Amperian_TRAP BC MIRRORED(probe)
function [I_amp, t_amp, tag, ProbeProfile, AmperProfile, ZerosProfile] = Calc_Amperian_TrapBCmirrored(probe)
% This version of Amperian assumes a linearly varying field between probes 
% for integration, using a trapezoidal rule, but also uses mirroring and 
% boundary conditions to account for regions where there is poor coverage
% by the surface probe array.  Assumptions are made where there are no
% probes.
%
% The Amperian calculation takes in the probe positions along the machine
% perimeter (arc lengths) and then calculations the distances between the
% LIVE probes to determine their Amperian lengths.  Then it's simply:
% I_tor = 1/u_o * sum_over_probes(B_p*L_amperian).
global OS
global PREFIX

OS = 'WIN';
PREFIX = 'T:\';

    %The Profile structures to record the profiles used to calculate the current
    ProbeProfile.name      = 'Probes';                                     % Legend entry for plotting the profile
    ProbeProfile.ptype     = 'xk';                                         % Line type for plotting
    ProbeProfile.t         = probe(1).t_tot';                              % Time basis (to make a movie)
    ProbeProfile.ArcLength = [];                                           % Arclength used for x-axis
    ProbeProfile.Bpol      = [];                                           % B_pol used for y-axis

    ZerosProfile.name      = 'Zeros';                                      % Legend entry for plotting the profile
    ZerosProfile.ptype     = 'ok';                                         % Line type for plotting
    ZerosProfile.t         = probe(1).t_tot';                              % Time basis (to make a movie)
    ZerosProfile.ArcLength = [];                                           % Arclength used for x-axis
    ZerosProfile.Bpol      = [];                                           % B_pol used for y-axis    
    
    AmperProfile.name      = 'Amperian';                                   % Legend entry for plotting the profile
    AmperProfile.ptype     = '-r';                                         % Line type for plotting
    AmperProfile.t         = probe(1).t_tot';                              % Time basis (to make a movie)
    AmperProfile.ArcLength = [];                                           % Arclength used for x-axis
    AmperProfile.Bpol      = [];                                           % B_pol used for y-axis
        
    if strcmpi(OS,'WIN')
        filename = [PREFIX 'A_Probe_Calibration\Codes- Analysis\Plasma Current\TetmeshArcLengths.txt'];
    elseif strcmpi(OS,'LINUX')
        filename = [PREFIX 'A_Probe_Calibration/Codes- Analysis/Plasma Current/TetmeshArcLengths.txt'];
    end
    
    [probe_name, probe_position] = textread(filename, '%s %f', 'delimiter', ',');
    for i = 1:length(probe_name)                                           % Turn the two vectors into a struct indexed by probe name
        position.(char(probe_name{i})) = probe_position(i);                % This is the position in terms of Arc Length (m), right now it's just used for the FIN position.
    end                                                                    %   
    for i = 1:length(probe)                                                % Go through each live probe in the array
%       probe(i).L = position.(probe(i).name(3:5));                        %   and assign the probe its arc length
%       probe(i).L is written by the probe_readin from Taylor
        ProbeProfile.ArcLength = [ProbeProfile.ArcLength, probe(i).L];     % This is a matrix, with the rows corresponding to different times and the columns to different ArcLengths
        ProbeProfile.Bpol = [ProbeProfile.Bpol, probe(i).b_tot'];          % 
    end                                                                    % 
    
    %         meters     Zero-Field Locations at Inside Corners, provided with Chris Hansen out to a billion sig figs
%   corner = [0.041,...             % Y Small nose cone
    corner = [0.400655945898837,... % Small Cone to Large Plate
              0.545685914134230,... % Large Plate to Large Cone
              0.869671598004181,... % Gap Y side
              0.9215,...            % Gap X side (0.929611598004599) (.921 by my assessment)
              1.253597281874924,... % Large Cone to Large Plate (1.253597281874924) (1.25 by my assessment)
              1.39];                % Large Plate to Small Cone (1.398627249150018) 
%              1.758];              % X Small nose cone
            
    ZerosProfile.ArcLength = corner;                                         % Arclength used for x-axis
    ZerosProfile.Bpol      = zeros([length(probe(1).t_tot),length(corner)]); % B_pol used for y-axis    
    % Do the Amperian math & for a time basis: I_enclosed = B*dl/u_o

    p = 0;                                                                 % p = probe counter, indicating index of current probe
    c = 0;                                                                 % c = corner counter, indicating index of current corner
    
    if corner(1) < probe(1).L                                              % If we're putting a BC on the nosecone, then...
        this.pc = 'c';                                                     % Initialize the first B-field to be the corner on the nose of the small cone
        this.L  = corner(1);                                               % Location of this B-value in terms of Arc Length
        this.B  = zeros(size(probe(1).b_tot'));                            % Value of the B field
        this.i  = 1;                                                       % Index for which probe or corner you're on

        next.pc = 'p';                                                     % Initialize the next B-field to be the first probe on the small cone
        next.L  = probe(1).L;                                              % Location of the B-value in terms of Arc Length
        next.B  = probe(1).b_tot';                                         % Value of the B field at the location of the next point
        next.i  = 1;                                                       % Index for which probe or corner you're on
        
    else                                                                   % If we're not putting a BC on the nosecone, then...
        this.pc = 'p';                                                     % Initialize the first B-field to be the first live probe on the small cone
        this.L  = probe(1).L;                                              % Location of the B-value in terms of Arc Length
        this.B  = probe(1).b_tot';                                         % Value of the B field at the location of the current point
        this.i  = 1;                                                       % Index for which probe or corner you're on
        
        next.pc = 'p';                                                     % Initialize the second B-field to be the second live probe on the small cone
        next.L  = probe(2).L;                                              % Location of this B-value in terms of Arc Length
        next.B  = probe(2).b_tot';                                         % Value of the B field
        next.i  = 2;                                                       % Index for which probe or corner you're on
    end

    I_amp = zeros(size(this.B));
    t_amp = probe(1).t_tot;
%     disp('New Array')
    while p <= length(probe) && c <= length(corner)                        % Keep going until I get to the last probe  
%         disp(['this.pc = ' this.pc])
%         disp(['this.L  = ' num2str(this.L)])
%         disp(['this.i  = ' num2str(this.i)])
%         disp(['next.pc = ' next.pc])
%         disp(['next.L  = ' num2str(next.L)])
%         disp(['next.i  = ' num2str(next.i)])
        if (this.pc == 'c') && (next.pc == 'p')                            % If we're heading out of a corner, calculate the trapezoid and reflect back
            L     = next.L - this.L;                                       % Calculate the base of the trapezoid
            B_avg = (this.B + next.B)/2;                                   % Calculate the top of the trapezoid (averaged points)
            I_amp = I_amp + (L * B_avg / (4*pi*10^-7));                    % Add the trapezoidal segment to Itor integral
            dx = next.L - this.L;                                          % Calculate dx for use in the slope (float)
            dy = next.B - this.B;                                          % Calculate dy for use in the slope (vector, length(B)x1)
            if (this.i == 1) && (probe(1).L < corner(1))                   % If it's the first corner w/a BC on the nosecone, reflect back to zero
                Lr     = this.L;                                           % Calculate the reflection w/dx & dy & dL
                Br_avg = (0 + (dy/dx)*Lr)/2;                               % Basically you're calculating the triangle from the zero boundary condition point
                AmperProfile.Bpol      = (dy/dx)*Lr;                       % Make the first data point
                AmperProfile.ArcLength = 0;                                %
            else                                                           % If it's a later corner, reflect halfway back to the previous corner
                Lr     = (this.L - corner(this.i-1))/2;                    % Calculate the reflection length
                Br_avg = (0 + (dy/dx)*Lr)/2;                               % Calculate the B-field using the slope 
                AmperProfile.Bpol      = [AmperProfile.Bpol,     (dy/dx)*Lr];  % Make the first data point
                AmperProfile.ArcLength = [AmperProfile.ArcLength, this.L-Lr];  %                
            end                                                            % 
            I_amp  = I_amp + (Lr * Br_avg / (4*pi*10^-7));                 % 
            AmperProfile.Bpol      = [AmperProfile.Bpol,      this.B];     % This is a matrix, with the rows corresponding to different times 
            AmperProfile.ArcLength = [AmperProfile.ArcLength, this.L];     %  and the columns to different ArcLengths  
            
        elseif (this.pc == 'p') && (next.pc == 'p')                        % If we're moving between probes, simply calculate the trapezoid
            if (this.i == 1) && (corner(1) > probe(1).L)                   % If it's the first probe and there is no BC on the nosecone
                dx = probe(1).L + (position.FIN - probe(end).L);           % then we'll have to interpolate to zero arc length
                dy = probe(1).b_tot' - probe(end).b_tot';                  % Calcule dx and dy for use in the slope
                Lr = probe(1).L;                                           % 
                Br_avg = (probe(1).b_tot' + (probe(1).b_tot' - (dy/dx)*Lr))/2;  % 
                I_amp  = I_amp + (Lr * Br_avg / (4*pi*10^-7));             % 
                AmperProfile.Bpol = [probe(1).b_tot' - (dy/dx)*Lr];        % Put in the first data point
                AmperProfile.ArcLength = [0];                              %  at zero arc length
            end
                
            L     = next.L - this.L;                                       % if it's a couple of probes somewhere in the middle
            B_avg = (this.B + next.B)/2;                                   %
            I_amp = I_amp + (L * B_avg / (4*pi*10^-7));                    %
            AmperProfile.Bpol      = [AmperProfile.Bpol,      this.B];     % This is a matrix, with the rows corresponding to different times
            AmperProfile.ArcLength = [AmperProfile.ArcLength, this.L];     %  and the columns to different ArcLengths
                
            if (next.i == length(probe)) && (corner(1) > probe(1).L)       % If it's the last probe and there is no BC on the nosecone
                dx = probe(1).L + (position.FIN - probe(end).L);           % then we'll have to interpolate to the final arc length
                dy = probe(1).b_tot' - probe(end).b_tot';                  % Calcule dx and dy for use in the slope
                Lr = position.FIN - probe(end).L;                          %  
                Br_avg = (probe(end).b_tot' + ((dy/dx)*Lr + probe(end).b_tot'))/2;   % 
                I_amp  = I_amp + (Lr * Br_avg / (4*pi*10^-7));             % 
                AmperProfile.Bpol      = [AmperProfile.Bpol,      next.B, probe(end).b_tot' + (dy/dx)*Lr]; % Put in the first data point
                AmperProfile.ArcLength = [AmperProfile.ArcLength, next.L,                   position.FIN]; %  at zero arc length
            end
            
        elseif (this.pc == 'p') && (next.pc == 'c')                        % If we're heading into a corner, calculate the trapezoid and reflect forwards
            L     = next.L - this.L;                                       %
            B_avg = (this.B + next.B)/2;                                   % 
            I_amp = I_amp + (L * B_avg / (4*pi*10^-7));                    % 
            AmperProfile.Bpol      = [AmperProfile.Bpol,      this.B];     % This is a matrix, with the rows corresponding to different times
            AmperProfile.ArcLength = [AmperProfile.ArcLength, this.L];     %  and the columns to different ArcLengths
            dx = next.L - this.L;                                          %
            dy = next.B - this.B;                                          % 
            if (next.i == length(corner)) && (corner(1) < probe(1).L)      % If it's the last corner, and there are nosecone BCs, reflect forward to position.FIN
                Lr     = position.FIN - next.L;                            % 
                Br_avg = (0 + (dy/(-dx))*Lr)/2;                            % 
                AmperProfile.Bpol      = [AmperProfile.Bpol,      zeros(size(dy)),  -(dy/dx)*Lr];   % Make the first data point
                AmperProfile.ArcLength = [AmperProfile.ArcLength, next.L           position.FIN];  % 
            else                                                           % If it's a previous corner, reflect forwards to halfway to the next corner
                Lr     = (corner(next.i+1) - corner(next.i))/2;            %
                Br_avg = (0 + (dy/(-dx))*Lr)/2;                            % 
                AmperProfile.Bpol      = [AmperProfile.Bpol,      zeros(size(dy)), -(dy/dx)*Lr];  % Make the first data point
                AmperProfile.ArcLength = [AmperProfile.ArcLength, next.L,            next.L+Lr];  % 
            end                                                            % 
            I_amp  = I_amp + (Lr * Br_avg / (4*pi*10^-7));                 % 
            
        elseif (this.pc == 'c') && (next.pc == 'c')
            %Do nothing, it's already taken care of
        end
                
        % Advancement Code...
        if this.pc == 'c'                                                  % If the outgoing point was a corner, 
            c = this.i;                                                    %   then store the corner number to a holding variable c
        else                                                               % But if the outgoing point was a probe,
            p = this.i;                                                    %   then store the probe number to a holding variable p
        end                                                                %
        
        this.pc = next.pc;                                                 % Push out the outgoing (this) point and replace
        this.L  = next.L;                                                  % it with the old next point.  This will allow us to
        this.B  = next.B;                                                  % grab the new next point...
        this.i  = next.i;                                                  % 
        
        if (this.pc == 'p') && (probe(1).L < corner(1)) && ...             % If we're on the last probe, and there are no nose-cone BCs
                (this.i == length(probe))                                  % 
            p = 400;                                                       %  then give "p" a huge number to step us out of the while loop
        elseif (this.pc == 'c' && (probe(1).L) > corner(1)) && ...         % If we're on the last corner, and there ARE nose-cone BCs
                (this.i == length(corner))                                 % 
            c = 400;                                                       %  then give "c" a huge number to step us out of the while loop
        else                                                               % Otherwise, we're not at the end, so carry on
            if this.pc == 'p'                                              % If we are currently on a probe then the next probe is (this.i+1), and the next corner is (c+1)
                if c == length(corner)
                    next.pc = 'p';                                         %  then it's the next point and written to the "next" variable 
                    next.L  = probe(this.i+1).L;                           % 
                    next.B  = probe(this.i+1).b_tot';                      % 
                    next.i  = this.i + 1;                                  % 
                elseif probe(this.i+1).L < corner(c+1)                     % And if the probe is nearer in arc length L than the probe
                    next.pc = 'p';                                         %  then it's the next point and written to the "next" variable 
                    next.L  = probe(this.i+1).L;                           % 
                    next.B  = probe(this.i+1).b_tot';                      % 
                    next.i  = this.i + 1;                                  % 
                else                                                       % Otherwise, if a corner is nearer in arc length, L, to "this" point
                    next.pc = 'c';                                         %   then the corner is taken as the "next" point.
                    next.L  = corner(c+1);                                 % 
                    next.B  = zeros(size(probe(1).b_tot'));                %  
                    next.i  = c+1;                                         % Give the corner it's proper corner number
                end                                                        % 
            elseif this.pc == 'c'                                          % If we are currently on a corner, then the next probe is (p+1) and the next corner is (this.i+1)
                if this.i == length(corner)
                    next.pc = 'p';                                         %
                    next.L  = probe(p+1).L;                                %
                    next.B  = probe(p+1).b_tot';                           %
                    next.i  = p+1;                                         %
                elseif probe(p+1).L < corner(this.i+1)                     % 
                    next.pc = 'p';                                         %
                    next.L  = probe(p+1).L;                                %
                    next.B  = probe(p+1).b_tot';                           %
                    next.i  = p+1;                                         %
                else                                                       % 
                    next.pc = 'c';                                         % 
                    next.L  = corner(this.i+1);                            % 
                    next.B  = zeros(size(probe(1).b_tot'));                % 
                    next.i  = this.i+1;                                    % 
                end
            end
        end
        
    end
    
    % RETURN: [I_amp, t_amp, tag, ProbeProfile, AmperProfile]
    tag = ['\I_tor_SPA' probe(1).name((end-2):end)];
end


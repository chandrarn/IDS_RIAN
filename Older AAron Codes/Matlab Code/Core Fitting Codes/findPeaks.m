function[chan_span, param] = findPeaks(data, peaks, param)
% 'data' is the full movie over the time of the shot
% 'time' is the time base in [ms] corresponding to the movie
% 'peaks' is an ('n_chan' x 3) matrix containing channel number, x, and y
% coordinates
%
% 'chan_span' is an ('n_chan' x 3) matrix containing channel number, left
% bound, and right bound for each channel (left and right are visual
% descriptions refering to lower and higher x value on the CCD).
%
% UPDATE June 14th, 2012:
% Adding 'width' factor so that the binning width can be narrowed
%
width = 2; % Full width to bin in pixels.  If set to zero, reverts to
           % original design which bins half distance to each neighboring
           % channel
hlfwdth = width / 2; % half the width for ease of notation

if isnan(peaks) % this code will do its best without any guidance
    % first, find time point with maximum brightness
    
    sumOfFrame = sum(sum(data, 3), 2); % sums over two data dimensions of 'data' array leaving
                                       % scalar sum of all counts for each
                                       % time.
    [dummy, n_bright] = max(sumOfFrame); % identify time index with most light
    
    % find location of approximate center in wavelenth space
    
    [dummy, center] = max(sum(data(n_bright, :, :), 3));
    
    % INCOMPLETE, COME BACK TO THIS LATER
    
else % blindly follow instructions from 'peaks' array
    
    [param.n_chan, dummy] = size(peaks); % extract number of channels
    chan_span = zeros(param.n_chan, 3); % preallocate aray
    chan_span(:, 1) = peaks(:, 1); % fill in channel numbers
    
    chan_sp_ave = median(diff(peaks(:, 2))); % finds the median channel spacing
                                             % (a 'mean' would be skewed by
                                             % non-consecutive channels)
    
    % first channel
    if chan_span(2, 1) == chan_span(1, 1) + 1 % the next channel is consecutive
        half = (peaks(2, 2) - peaks(1, 2)) / 2; % half the distance to the next channel
        if or(~width, hlfwdth > half) % width=0 or the preset width is greater than 
                                      % half the dist. to next channel
            chan_span(1, 2) = peaks(1, 2) - half; % lower (left) bound on channel domain
            chan_span(1, 3) = peaks(1, 2) + half; % upper (right) bound on channel domain
        else % width is specified nonzero and smaller than half the dist. to next channel
            chan_span(1, 2) = peaks(1, 2) - hlfwdth; % lower (left) bound on channel domain
            chan_span(1, 3) = peaks(1, 2) + hlfwdth; % upper (right) bound on channel domain
        end
    else % next channel is NOT consecutive
        if ~width % width is not specified
            chan_span(1, 2) = peaks(1, 2) - (chan_sp_ave / 2); % channel domain based solely on median
            chan_span(1, 3) = peaks(1, 2) + (chan_sp_ave / 2);
        else
            chan_span(1, 2) = peaks(1, 2) - hlfwdth; % channel domain based solely on specified width
            chan_span(1, 3) = peaks(1, 2) + hlfwdth;
        end
    end
    
    % loop through middle channels
    for n = 2:(param.n_chan - 1)
        if and(chan_span(n, 1) ~= chan_span(n-1, 1) + 1, chan_span(n, 1) ~= chan_span(n+1, 1) - 1)
            % the channel is totally isolated with no consecutive neighbors
            if ~width % width is not specified
                chan_span(n, 2) = peaks(n, 2) - (chan_sp_ave / 2); % channel domain based solely on median
                chan_span(n, 3) = peaks(n, 2) + (chan_sp_ave / 2);
            else
                chan_span(n, 2) = peaks(n, 2) - hlfwdth; % channel domain based solely on specified width
                chan_span(n, 3) = peaks(n, 2) + hlfwdth;
            end
        elseif and(chan_span(n, 1) == chan_span(n-1, 1) + 1, chan_span(n, 1) ~= chan_span(n+1, 1) - 1)
            % the channel has a left neighbor but NOT a right neighbor
            half = (peaks(n, 2) - peaks(n-1, 2)) / 2; % half the distance to the left neighbor
            if or(~width, hlfwdth > half) % width=0 or the preset width is greater than 
                                          % half the dist. to left neighbor channel
                chan_span(n, 2) = peaks(n, 2) - half;
                chan_span(n, 3) = peaks(n, 2) + half;
            else % width is specified nonzero and smaller than half the dist. to left neighbor channel
                chan_span(n, 2) = peaks(n, 2) - hlfwdth; % channel domain based on specified width
                chan_span(n, 3) = peaks(n, 2) + hlfwdth;
            end
        elseif and(chan_span(n, 1) ~= chan_span(n-1, 1) + 1, chan_span(n, 1) == chan_span(n+1, 1) - 1)
            % the channel does NOT have a left neighbor, but HAS a right
            % neighbor
            half = (peaks(n+1, 2) - peaks(n, 2)) / 2; % half the distance to the right neighbor
            if or(~width, hlfwdth > half) % width=0 or the preset width is greater than 
                                          % half the dist. to right neighbor channel
                chan_span(n, 2) = peaks(n, 2) - half;
                chan_span(n, 3) = peaks(n, 2) + half;
            else % width is specified nonzero and smaller than half the dist. to right neighbor channel
                chan_span(n, 2) = peaks(n, 2) - hlfwdth; % channel domain based on specified width
                chan_span(n, 3) = peaks(n, 2) + hlfwdth;
            end
        else % the channel must have both neighbors
            if or(~width, or(hlfwdth > (peaks(n, 2) - peaks(n-1, 2)) / 2, hlfwdth > (peaks(n+1, 2) - peaks(n, 2)) / 2))
                % width is not specified or specified width is greater than
                % half distance to either neighbor
                chan_span(n, 2) = peaks(n, 2) - (peaks(n, 2) - peaks(n-1, 2)) / 2; 
                                    % half the distance to the left neighbor
                chan_span(n, 3) = peaks(n, 2) + (peaks(n+1, 2) - peaks(n, 2)) / 2; 
                                    % half the distance to the right
                                    % neighbor
            else % width is specified nonzero and smaller than half the dist. either neighbor
                chan_span(n, 2) = peaks(n, 2) - hlfwdth; % channel domain based on specified width
                chan_span(n, 3) = peaks(n, 2) + hlfwdth;
            end
        end
    end
    
    % last channel
    if chan_span(param.n_chan, 1) == chan_span(param.n_chan-1, 1) + 1 % the last channel has a left neighbor
        half = (peaks(param.n_chan, 2) - peaks(param.n_chan-1, 2)) / 2; % half the distance to the previous channel
        if or(~width, hlfwdth > half) % width=0 or the preset width is greater than 
                                      % half the dist. to left neighbor channel
            chan_span(param.n_chan, 2) = peaks(param.n_chan, 2) - half; % lower (left) bound on channel domain
            chan_span(param.n_chan, 3) = peaks(param.n_chan, 2) + half; % upper (right) bound on channel domain
        else % width is specified nonzero and smaller than half the dist. to left neighbor channel
            chan_span(param.n_chan, 2) = peaks(param.n_chan, 2) - hlfwdth; % lower (left) bound on channel domain
            chan_span(param.n_chan, 3) = peaks(param.n_chan, 2) + hlfwdth; % upper (right) bound on channel domain
        end
    else % the last channel is isolated
        if ~width % width is not specified
            chan_span(param.n_chan, 2) = peaks(param.n_chan, 2) - (chan_sp_ave / 2); % channel domain based solely on median
            chan_span(param.n_chan, 3) = peaks(param.n_chan, 2) + (chan_sp_ave / 2);
        else % width is specified
            chan_span(param.n_chan, 2) = peaks(param.n_chan, 2) - hlfwdth; % channel domain based on 
            chan_span(param.n_chan, 3) = peaks(param.n_chan, 2) + hlfwdth;
        end
    end
        
end
end
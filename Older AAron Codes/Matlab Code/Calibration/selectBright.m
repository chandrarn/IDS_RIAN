function[data3, time] = selectBright(data2, time)
% Parameters
lamp_freq = 58; % Hz
max_evals = 200; % maximum evaluations for fminsearch

% first, make single vector with all light.

sumOfFrame = sum(sum(data2, 3), 2); % sums over two data dimensions of 'data' array leaving
                                   % scalar sum of all counts for each
                                   % time.
                                   
%% Fit function to emission

% Parameter first guesses for rectified sine fit

A = zeros(1, 3); % holder for fit parameters
A(1) = max(sumOfFrame); % Amplitude in [counts]
A(2) = lamp_freq * 2 * pi / 1000; % frequency in [rad/ms]
A(3) = A(2)*time(find(sumOfFrame > 0.9*A(1), 1)); % Phase offset

% Find best fit

[estimates] = optimize(time, sumOfFrame, A, max_evals); % call fit function
A = estimates;
% fit = A(1) .* abs(sin(A(2) .* time + A(3))); % evaluate best fit.

% Find ideal times of maximum light
tmin = (pi/2 - A(3))/A(2);
T = 2*pi/A(2); % period
n_pos = floor((time(end) - tmin)/T); % number of positive peak emission times
t_peaks = tmin + T*[1:n_pos]; % times of ideal peaks

% loop to identify closest frames to peak light
n_peaks = [];
for n = 1:length(t_peaks)
    [dummy, Ind] = min(abs(t_peaks(n) - time));
    n_peaks = [n_peaks, Ind]; % data time points of maxmum light
end

% Select only brightest data
data3 = data2(n_peaks, :, :);
time = time(n_peaks);

end

function[estimates] = optimize(time, y, A, max_evals)
options = optimset('MaxFunEvals', max_evals);
model = @crazy_func;
estimates = fminsearch(model, A);

    function [sse, y_fit] = crazy_func(A)
        y_fit = A(1) .* abs(sin(A(2) .* time + A(3)));
        sse = sum((y_fit - y).^2);
    end
end
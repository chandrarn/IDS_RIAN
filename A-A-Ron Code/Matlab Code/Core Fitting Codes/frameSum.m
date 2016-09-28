function[time, data2] = frameSum(time, data2, frame_sum, param)

[n_time, ymax, xmax] = size(data2); % determine size of data array

n_time = floor(n_time / frame_sum); % resulting number of time points
data3 = zeros(n_time, ymax, param.n_chan); % preallocate new data array

% Interpolate time base

dt = time(2) - time(1); % original time step
time = time(1:frame_sum:frame_sum * n_time) + (dt * (frame_sum - 1) / 2);

% Sum data

for n = 1:n_time % loop over all new time points
    data3(n, :, :) = sum(data2(n * frame_sum - frame_sum + 1:n * frame_sum, :, :), 1);
end

data2 = data3;
end
function[temp, vel] = sanitizeOut(temp, vel, limits)

[n_time, n_chan] = size(temp);

for n = 1:n_time
    for m = 1:n_chan
        if or(temp(n, m) > limits(1, 2), temp(n, m) < limits(1, 1))
            temp(n, m) = NaN;
        end
        if or(vel(n, m) > limits(2, 2), vel(n, m) < limits(2, 1))
            vel(n, m) = NaN;
        end
    end
end
end
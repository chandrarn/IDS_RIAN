function dat = trimRange(dat, chan_range)

% trims all parameters in 'dat' structure to conform to channel range

k = 1;
for n = 1:length(chan_range)
    try
        ind(k) = find(dat.peaks == chan_range(n));
        k = k+1;
    end
end

dat.impacts = dat.impacts(ind);
dat.temp = dat.temp(:, ind);
dat.vel = dat.vel(:, ind);
dat.peaks = dat.peaks(ind);
try % these won't exist for NIMROD results
    dat.int = dat.int(:, ind);
    dat.fit_par = dat.fit_par(:, ind, :);
    dat.bounds = dat.bounds(:, ind, :);
    dat.guesses = dat.guesses(:, ind, :);
end
end
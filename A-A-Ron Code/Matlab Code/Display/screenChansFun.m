function[dat] = screenChansFun(dat, scChans)

[n_time, n_chan] = size(dat.temp);

a = length(scChans);
wipe = NaN * zeros(n_time, 1);

for n = 1:a
    b = find(dat.peaks == scChans(n));
    if b % my desired channel DOES exist
        dat.temp(:, b) = wipe;
        dat.vel(:, b) = wipe;
    end
end

end
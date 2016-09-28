function[dat] = screenAreaFun(dat, scArea)

wipe = find(squeeze(dat.fit_par(:, 1, :)) < scArea);

dat.temp(wipe(:)) = NaN * ones(length(wipe(:)), 1);
dat.vel(wipe(:)) = NaN * ones(length(wipe(:)), 1);

end
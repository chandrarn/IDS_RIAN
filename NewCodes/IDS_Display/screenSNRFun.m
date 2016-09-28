function[dat] = screenSNRFun(dat, scSNR)

wipe = find(dat.residual ./ squeeze(dat.fit_par(:, 1, :)) > scSNR);

dat.temp(wipe(:)) = NaN * ones(length(wipe(:)), 1);
dat.vel(wipe(:)) = NaN * ones(length(wipe(:)), 1);

end
function[dat] = screenIntFun(dat, scInt)

wipe = find(dat.int(:) < scInt);

dat.temp(wipe(:)) = NaN * ones(length(wipe(:)), 1);
dat.vel(wipe(:)) = NaN * ones(length(wipe(:)), 1);

end
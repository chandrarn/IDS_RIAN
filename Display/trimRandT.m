function dat = trimRandT(dat, chan_range,time_range)
% Rian Chandra: Updated to trim time as well.
% trims all parameters in 'dat' structure to conform to channel range

if ~length(time_range) % if only trimming channels, use all time
    tNorm = [1:length(dat.time)];
    tInj = [1:length(dat.iinjxTime)];
    tTor = [1:length(dat.ItorTime)];
else
    tNorm = timeRange(time_range,dat.time); % normal time
    tInj = timeRange(time_range,dat.iinjxTime); % injector time
    tTor = timeRange(time_range,dat.ItorTime); % torroidal time
end

k = 1;
for n = 1:length(chan_range)
    try
        ind(k) = find(dat.peaks == chan_range(n));
        k = k+1;
    end
end
size(dat.impacts);
dat.impacts = dat.impacts(ind);
dat.temp = dat.temp(tNorm, ind);
dat.vel = dat.vel(tNorm, ind);
dat.peaks = dat.peaks(ind);
try % these won't exist for NIMROD results
    dat.int = dat.int(tNorm, ind);
    dat.fit_par = dat.fit_par(tNorm, ind, :);
    dat.bounds = dat.bounds(tNorm, ind, :);
    dat.guesses = dat.guesses(tNorm, ind, :);
    dat.time = dat.time(tNorm);
    dat.iinjxTime = dat.iinjxTime(tInj);
    dat.iinjx = dat.iinjx(tInj);
    dat.iinjyTime = dat.iinjyTime(tInj);
    dat.iinjy = dat.iinjy(tInj);
    dat.ItorTime = dat.ItorTime(tTor);
    dat.Itor = dat.Itor(tTor);
    dat.tempU = dat.tempU(tNorm, ind, :); % Need to go at bottom: in cases
    dat.tempL = dat.tempL(tNorm, ind, :); % where batch correct did not produce
    dat.velU = dat.velU(tNorm, ind, :);   % error files, the try will exit, and
    dat.velL = dat.velL(tNorm, ind, :);   % the rest of the try will not execute.
    dat.dparam = dat.dparam(tNorm, ind, :);
    dat.stddev = dat.stddev(tNorm, ind, :);
    
end
end


function tOut = timeRange(tIn,tBase)
% find time bounds
xMin = find(tBase>=tIn(1));
xMin = xMin(1);
xMax = find(tBase>=tIn(2));
xMax = xMax(1);

tOut = xMin:xMax;
end

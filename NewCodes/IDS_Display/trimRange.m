function dat = trimRange(dat, chan_range, error, time_range, trimParam)

% trims all parameters in 'dat' structure to conform to channel range

% Find time limits, in units of dat.time
if isempty(time_range)
    indT = 1:length(dat(1).time);
else
    [~, n1] = min(abs(dat(1).time - time_range(1)));
    [~, n2] = min(abs(dat(1).time - time_range(2)));
    indT = n1:n2;
end

k = 1;
for n = 1:length(chan_range)
    try
        ind(k) = find(dat(1).peaks == chan_range(n));
        k = k+1;
    end
end

if (trimParam)
    dat(1).param.peaks = dat(1).param.peaks(ind,:);
    dat(1).param.PIX_SP = dat(1).param.PIX_SP(ind);
    dat(1).param.REL_INT = dat(1).param.REL_INT(ind);
    dat(1).param.impacts = dat(1).param.impacts(ind);
    dat(1).param.Center = dat(1).param.Center(ind,:);
    dat(1).param.Inst_Temp = dat(1).param.Inst_Temp(ind,:);
end
dat(1).peaks = dat(1).peaks(ind);
dat(1).impacts = dat(1).impacts(ind);

dat(1).time = dat(1).time(indT);
for n = 1:length(dat)
    dat(n).temp = dat(n).temp(indT, ind);
    dat(n).vel = dat(n).vel(indT, ind);
    dat(n).int = dat(n).int(indT, ind);
    dat(n).fit_par = dat(n).fit_par(indT, ind, :);
    dat(n).bounds = dat(n).bounds(indT, ind, :);
    dat(n).guesses = dat(n).guesses(indT, ind, :);
    if error
        dat(n).intL = dat(n).intL(indT, ind);
        dat(n).intU = dat(n).intU(indT, ind);
        dat(n).velL = dat(n).velL(indT, ind);
        dat(n).velU = dat(n).velU(indT, ind);
        dat(n).tempL = dat(n).tempL(indT, ind);
        dat(n).tempU = dat(n).tempU(indT, ind);
    end
end
end
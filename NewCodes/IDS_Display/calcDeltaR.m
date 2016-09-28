function dat = calcDeltaR(dat, uniformTemp, t_avg)

temp1 = dat.temp - uniformTemp; % subtract off a uniform temperature

for n = 1:length(temp1(:))
    if temp1(n) < 0 % call any negative temperature zero
        temp1(n) = 0;
    end
end

v = sqrt(1.6e-19 .* temp1 ./ dat.param.IonMass);

dat.deltaR = 1e2 * 1e-6 * t_avg * v; % convert t_avg from [us] to [s]
% and [m] to [cm]

end
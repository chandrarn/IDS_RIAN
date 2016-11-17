function y = singletGauss1DLM(x, par, c)
%
% par(1) = area (the function is a normal distribution)
% par(2) = x0
% par(3) = sigx
% par(4) = offset
%
% x is a column vector

y = par(1)/(sqrt(2*pi)*par(3)) * exp(-0.5 * (x-par(2)).^2/(par(3).^2)) + par(4);
end
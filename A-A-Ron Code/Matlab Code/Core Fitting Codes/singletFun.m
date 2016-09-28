function y = singletFun(par, xdata)
%
% par(1) = area (normal distribution)
% par(2) = y0 (wavelength initial center guess)
% par(3) = sigx
% par(4) = offset
%
% xdata is a vector of points to fit the data onto.
y = par(1) / (par(3) * 2.5066) * exp(-(xdata - par(2)).^2 / (2*par(3)^2)) ...
    + par(4);

% NB: '2.5066' is actually 'sqrt(2*pi)' evaluated to make code faster
end
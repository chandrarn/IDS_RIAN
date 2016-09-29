function z = singletGauss2DLM(x, par, c)
%
% par(1) = "volume" (the function is a normal distribution)
% par(2) = x0
% par(3) = y0
% par(4) = sigx
% par(5) = sigy
% par(6) = offset
%
% x is a 2 column array where x(:, 1) is mesh X and x(:, 2) is mesh Y

z = par(1)/(2*pi*par(4)*par(5)) * exp(-0.5 * (((x(:,1) - par(2)) ./ par(4)).^2 + ...
    ((x(:,2) - par(3)) ./ par(5)).^2)) + par(6);
%z = z';
end
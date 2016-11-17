function y = SineFitLM(x, par, c)
%
% par(1) = "volume" (the function is a normal distribution)
% par(2) = x0
% par(3) = y0
% par(4) = sigx
% par(5) = sigy
% par(6) = offset
%
% x is a 2 column array where x(:, 1) is mesh X and x(:, 2) is mesh Y

y = par(1) + par(2)*sin(par(3)+x*par(4)*2*pi);
end
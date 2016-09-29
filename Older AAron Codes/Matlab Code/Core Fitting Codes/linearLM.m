function z = linearLM(x, par, c)
%
% par(1) = slope
% par(2) = offset
%
% x is a 1 column vector to fit over

z = par(1).*x +par(2);
%z = z';
end
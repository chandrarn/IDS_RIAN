function z = fitSine(par,x)

z = par(1)*sin(par(2)*pi*x + par(3)*pi) + par(4);
end
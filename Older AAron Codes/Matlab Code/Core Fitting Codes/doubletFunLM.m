function y = doubletFunLM(xdata, guess, c)
% c(1) = param.deltaPix(m)
% c(2) = param.ratio

y = guess(1) / (guess(3) * 2.5066) * (exp(-(xdata - guess(2)).^2 / (2*guess(3)^2)) + ...
    c(2) * exp(-(xdata - (guess(2) + c(1))).^2 / (2*guess(3)^2))) ...
    + guess(4);
y = y';

% NB: '2.5066' is actually 'sqrt(2*pi)' evaluated to make code faster
end
function y = doubletFun(guess, xdata, deltaPix, ratio)

y = guess(1) / (guess(3) * 2.5066) * (exp(-(xdata - guess(2)).^2 / (2*guess(3)^2)) + ...
    ratio * exp(-(xdata - (guess(2) + deltaPix)).^2 / (2*guess(3)^2))) ...
    + guess(4);

% NB: '2.5066' is actually 'sqrt(2*pi)' evaluated to make code faster
end
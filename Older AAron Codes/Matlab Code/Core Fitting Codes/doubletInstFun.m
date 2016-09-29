function y = doubletInstFun(guess, xfine, param, ch)

denom = param.FWHM(ch)^2 / (4*log(2)*param.PIX_SP(ch)^2); % denominator of Gaussian exponent

y = guess(1) / (guess(3) * 2.5066) * (exp(-(xfine - param.Center(ch)).^2 / denom) + ...
    param.ratio * exp(-(xfine - (param.Center(ch) + param.deltaPix(ch))).^2 / denom)) ...
    + guess(4);
end
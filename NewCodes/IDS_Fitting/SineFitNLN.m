function y = SineFitNLN(par,x, c)
%
% par(1) = Offset
% par(2) = Amplitude
% par(3) = Phase
% par(4) = Frequency
%
% x is the time vector, y is the sine wave
y = par(1) + par(2)*sin(par(3)+x*par(4)*2*pi);
end
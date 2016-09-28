%% [signal, t] = Preliminary_Filter(signal, t, n)
function [signal, t] = Preliminary_Filter(signal, t, n)
% Preliminary Filter is used to remove DC and slant offsets to improve 
% signals for subsequent integration.  It fits a line:
% written by JSW
%
% n = 0, y=b         horizontal line
% n = 1, y=mx+b      straigth line
% n = 2, y=ax^2+bx+c parabolic
% 
% Butterworth filtering didn't work well, and Read_In takes care of 
% baseline subtraction.  What we've settled on is 
% Poly fitting (using y=mx+b) to the pre and post plasma signal sections, 
% and subtracting off that fitted line.
% Removes any DC offset
% Removes any linear drift

global OUTPUT;                                                             %

[value,start] = min(abs(t- -2e-4));                                    % Find the pre-shot section
[value,stop]  = min(abs(t- 15e-3));                                    % Find the post-shot section
offset_coeff  = polyfit([t(1:start);t(stop:end)],[signal(1:start);signal(stop:end)],n);
signal_offset = signal-polyval(offset_coeff,t);                        %

signal = signal_offset;
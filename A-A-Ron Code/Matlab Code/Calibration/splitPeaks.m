function[peaks] = splitPeaks(peaks, isUpper)
% Chops off 'peaks' to match fibers being used for particular calibration
lastUpper = 36; % channel number of last fiber in upper array

if isUpper
    peaks = peaks(find(peaks(:, 1) <= lastUpper), :);
else
    peaks = peaks(find(peaks(:, 1) > lastUpper), :);
end
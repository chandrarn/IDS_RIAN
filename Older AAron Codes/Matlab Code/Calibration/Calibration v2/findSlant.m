function slope = findSlant(PEAKS)

% Finds 'm*x + b' slope parameter for adjusting centering estimates in motor
% calibration.  'x' is the wavelength space coordinate.

% Find inverse slope - if line is slanted slightly down to the right on the 
% image, this gives a very large positive slope.

m_1 = (PEAKS(end, 2) - PEAKS(1, 2)) / (PEAKS(1, 3) - PEAKS(end, 3));

slope = 1 / m_1; % 

end
% find the phase based on the param output from sine_fit, centers the
% phasing around 1.7 ms (easier to visually verify veracity)

% References off of the zero crossing to the LEFT of the 1.7 ms line

function phase = find_Phase(param)

timePoint = 0;.0017; % reference timepoint. Note, for 14500Hz, the closest true
% zero crossing is .00172414.

 initPhase = asin(sin(2*pi*(.0017)*param(5) + param(4))); % initial phase, -pi/2 - pi/2
 
 % is the sign of the derivative is positive, we're in the range where asin
 % is valid
 % if its zero (somehow) asin will stil be correct in reporting +- pi/2
 if sign(cos(2*pi*(.0017)*param(5) + param(4))) >=0 
     phase = initPhase;
     return;
 
 % if sign of derivative is negative, 
 elseif sign(cos(2*pi*(.0017)*param(5) + param(4))) == -1
     s = sign(initPhase);
     phase = s*(pi-abs(initPhase));
    return;
 end
end
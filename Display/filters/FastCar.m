% FastCar(data, pointsBefore, pointsAfter) is the boxcar average of data,
% using the passed number of points before and after.
%
% If only two arguments are passed, it splits it into a symmetric boxcar.
%
% The use of two arguments allows for even-width, non-prime width, non-symmetric
% boxcars, while with one argument this is prevented.
function smoothed = FastCar(data, pointsBefore, pointsAfter)
    % Argument Handling:
    if nargin < 3
        if mod(pointsBefore, 2) ~= 1
            error('FastCar:evenBoxcarWidth',...
                'You tried to make an even-width boxcar.  Stop that!');
        end
        if ~isprime(pointsBefore)
            error('FastCar:nonPrimeWidth',...
                'You tried to make your boxcar have a non-prime width.  You are clearly not a jedi, yet.');
        end
        pointsAfter = (pointsBefore - 1) / 2;
        pointsBefore = pointsAfter;
    end
    
    % Actual smoothing:
    smoothed = data;
    smoothed(pointsBefore + 1) = sum(data(1 : pointsBefore + 1 + pointsAfter)) / (pointsBefore + 1 + pointsAfter);
    for index = pointsBefore + 2 : length(data) - pointsAfter
        smoothed(index) = smoothed(index - 1) + (data(index + pointsAfter) - data(index - pointsBefore - 1)) / (pointsBefore + 1 + pointsAfter);
    end
end
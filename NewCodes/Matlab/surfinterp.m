function handle = surfinterp(varargin)
% This interpolates and plots surf data. I'm writing this because 'shading
% interp' is broken/ignored on my machine. Therefore a linear interpolation
% will do the same thing...
if nargin == 1
    handle = surf(varargin{1}');
    return;
elseif nargin == 2
    Z = varargin{1};
    factor = varargin{2};
    [X, Y] = ndgrid(1:size(Z, 1), 1:size(Z, 2));
elseif nargin == 3
    handle = surf(varargin{1}', varargin{2}', varargin{3}');
    return;
elseif nargin == 4
    X = varargin{1};
    Y = varargin{2};
    Z = varargin{3};
    factor = varargin{4};
else
    warning('too many input arguments for function ''surfinterp''')
end
    

F = griddedInterpolant(X, Y, Z);

x2 = linspace(X(1, 1), X(end, 1), factor * size(X, 1));
y2 = linspace(Y(1, 1), Y(1, end), factor * size(Y, 2));

[X2, Y2] = ndgrid(x2, y2);

handle = surf(X2', Y2', F(X2, Y2)');
function[Xf, Yf, Zf] = reconstructCalFits(par, brightWing)

xmin = round(par(1, 2)) - 2; % minimum x value of fit
xmax = round(par(end, 2)) + 2; % maximum x value of fit

ymin = round(mean(par(:, 3))) - brightWing; % minimum y value of fit
ymax = round(mean(par(:, 3))) + brightWing; % maximum y value of fit

[Xf, Yf] = meshgrid(xmin:0.1:xmax, ymin:0.1:ymax);

% reshape for function execution

xf(:, 1) = Xf(:);
xf(:, 2) = Yf(:);
zf = 0 * xf(:, 1); % initialize data vector as zeros

for n = 1:size(par, 1) % loop over all channels
    
    fit = singletGauss2D(par(n, :), xf);
    
    zf = zf + fit; % add new fit 
    
end

Zf = reshape(zf, size(Xf, 1), size(Xf, 2)); % reshape into 2D image
end
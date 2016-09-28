%% TEST OF LM-METHOD FITTING IN 2-D
function [p_fit,Chi_sq,sigma_p,sigma_y,corr,R2,cvg_hst] = Lm2DTest
close all; 
clear all;
shot = 12949910;% shot number
eval(sprintf('load(''dat%0.0f'');', shot)); % get dat

n=196;%timestamp
m=14;%chan 12

dp = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001];

x0 = round(dat.param.peaks(m, 2));
y0 = round(dat.param.Center(m));
        
xBound = x0 - dat.param.xWing : x0 + dat.param.xWing;
yBound = y0 - dat.param.yWing : y0 + dat.param.yWing;
disp('SHALALALALALALALALALALALA');

%[X, Y] = meshgrid(xBound, yBound);
%Z = squeeze(data(n, yBound, xBound));


[X,Y]=meshgrid(-5:.3:5,-5:.3:5);


% 
% clear all;
% [X, Y] = meshgrid(-5:.1:5, -5:.1:5);

 x(:, 1) = X(:);
 x(:, 2) = Y(:);
[sx sy]=size(X);
 Z = lm_func2d(x,[.1,1.7],0); % create some data
 %p_true = [1243.04667365151,97.8073597485189,40.2152933825534,1.50572796181989,2.14437493106966,28.2643984590859];
 %Z = singletGauss2DLM(x,p_true,0); % create some data
 

 x(:,3) = Z(:);
  figure; surf(X,Y,reshape(real(Z),sx,sy)); shading interp; hold on;
 t = [x(:,1) x(:,2) ];
%p_min =  [1100.04667365151,90.8073597485189,30.2152933825534,1.00572796181989,1.54437493106966,20.2643984590859]; % minimum expected parameter values
%p_max = [1300.04667365151,105.8073597485189,50.2152933825534,2.0572796181989,2.74437493106966,35.2643984590859]; % maximum expected parameter values
%p_init = [1200.04667365151,97.8073597485189,40.2152933825534,1.10572796181989,2.64437493106966,28.2643984590859]; % initial guess for parameter values
z_dat= x(:,3);
%Z= singletGauss2DLM(x,p_init,0);

%plot3(X,Y,reshape(Z,sx,sy),'*');

% [p_fit,Chi_sq,sigma_p,sigma_y,corr,R2,cvg_hst] = ...
% lm(@singletGauss2DLM, p_init, t, z_dat, .001,.001,p_min,p_max,0)
[p_fit,Chi_sq,sigma_p,sigma_y,corr,R2,cvg_hst] = ...
lm(@lm_func2d, [.5,1.1], t, z_dat, 1,.01,[0 0],[2 2],0)
%lm(func, p, t, y_dat, weight, dp, p_min, p_max, c)

Z=lm_func2d(t,p_fit,0);
plot3(X,Y,reshape(Z,sx,sy),'*');


end

function z_hat = lm_func2d(t,p,c)
    % example function used for nonlinear least squares curve-fitting
    % to demonstrate the Levenberg-Marquardt function, lm.m,
    % in two fitting dimensions
    x_dat = t(:,1);
    y_dat = t(:,2);
    z_hat = real(( p(1)*x_dat.^p(2) + (1-p(1))*y_dat.^p(2) ).^(1/p(2)));
end


function z = singletGauss2DLM(x, par, c)
%
% par(1) = "volume" (the function is a normal distribution)
% par(2) = x0
% par(3) = y0
% par(4) = sigx
% par(5) = sigy
% par(6) = offset
%
% x is a 2 column array where x(:, 1) is mesh X and x(:, 2) is mesh Y
%par(1)/(2*pi*par(4)*par(5)) * 
%1718.78480554029/(2*pi*1.97978544738108*2.11506601933119)
z =  par(1)/(2*pi*par(4)*par(5)) *  exp(-0.5 * ((((x(:,1)) - par(2)) ./ par(4)).^2 + (((x(:,2)) - par(3)) ./ par(5)).^2)) + par(6);
%z = z';
end

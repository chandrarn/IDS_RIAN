
function [alpha,beta,Chi_sq,y_hat,dydp] = lm_matx(func,t,p,y_dat,weight_sq,dp,c)
% [alpha,beta,Chi_sq,y_hat,dydp] = lm_matx(func,t,p,y_dat,weight_sq,{da},{c})
%
% Evaluate the linearized fitting matrix, alpha, and vector beta, 
% and calculate the Chi-squared error function, Chi_sq 
% Used by Levenberg-Marquard algorithm, lm.m   
% -------- INPUT VARIABLES ---------
% func  = function ofpn independent variables, p, and m parameters, p,
%         returning the simulated model: y_hat = func(t,p,c)
% t     = m-vectors or matrix of independent variables (used as arg to func)
% p     = n-vector of current parameter values
% y_dat = n-vector of data to be fit by func(t,p,c)  
% weight_sq = square of the weighting vector for least squares fit ...
%	    inverse of the standard measurement errors
% dp = fractional increment of 'p' for numerical derivatives
%      dp(j)>0 central differences calculated
%      dp(j)<0 one sided differences calculated
%      dp(j)=0 sets corresponding partials to zero; i.e. holds p(j) fixed
%      Default:  0.001;
% c  = optional vector of constants passed to y_hat = func(t,p,c)
%---------- OUTPUT VARIABLES -------
% alpha	= linearized Hessian matrix (inverse of covariance matrix)
% beta  = linearized fitting vector
% Chi_sq = 2*Chi squared criteria: weighted sum of the squared residuals WSSR
% y_hat = model evaluated with parameters 'p'

%   Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. November 2005
%   modified from: ftp://fly.cnuce.cnr.it/pub/software/octave/leasqr/
%   Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.

 global func_calls

 Npnt = length(y_dat);		% number of data points
 Npar = length(p);		% number of parameters 

 if nargin < 6
      dp = 0.001;
 end

 alpha = zeros(Npar);
 beta  = zeros(Npar,1);

 y_hat = feval(func,t,p,c);	% evaluate model using parameters 'p'
 func_calls = func_calls + 1;
size(y_dat);
size(y_hat);
 delta_y = y_dat - y_hat;	% residual error between model and data

 dydp = lm_dydp(func,t,p,y_hat,dp,c);

 alpha = dydp' * ( dydp .* ( weight_sq * ones(1,Npar) ) );  

 beta  = dydp' * ( weight_sq .* delta_y );
 
 Chi_sq = delta_y' * ( delta_y .* weight_sq ); 	% Chi-squared error criteria

% endfunction  # ------------------------------------------------------ LM_MATX


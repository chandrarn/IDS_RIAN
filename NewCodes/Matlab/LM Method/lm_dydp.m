
function dydp = lm_dydp(func,t,p,y,dp,c)
% dydp = lm_dydp(func,t,p,y,{dp},{c})
%
% Numerical partial derivatives (Jacobian) dy/dp for use with lm.m
% Requires n or 2n function evaluations, n = number of nonzero values of dp
% -------- INPUT VARIABLES ---------
% func = function of independent variables, 't', and parameters, 'p',
%        returning the simulated model: y_hat = func(t,p,c)
% t  = m-vector of independent variables (used as arg to func) 
% p  = n-vector of current parameter values
% y  = func(t,p,c) n-vector initialised by user before each call to lm_dydp
% dp = fractional increment of p for numerical derivatives
%      dp(j)>0 central differences calculated
%      dp(j)<0 one sided differences calculated
%      dp(j)=0 sets corresponding partials to zero; i.e. holds p(j) fixed
%      Default:  0.001;
% c  = optional vector of constants passed to y_hat = func(t,p,c)
%---------- OUTPUT VARIABLES -------
% dydp = Jacobian Matrix dydp(i,j)=dy(i)/dp(j)	i=1:n; j=1:m 

%   Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. November 2005
%   modified from: ftp://fly.cnuce.cnr.it/pub/software/octave/leasqr/
%   Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.


 global func_calls

 m=length(y);			% number of data points
 n=length(p);			% number of parameters

 if nargin < 5
	dp = 0.001*ones(1,n);
 end

 ps=p; dydp=zeros(m,n); del=zeros(n,1);         % initialize Jacobian to Zero

 for j=1:n                       % loop over all parameters

     del(j) = dp(j) * (1+abs(p(j)));  % parameter perturbation
     p(j)   = ps(j) + del(j);	      % perturb parameter p(j)

     if del(j) ~= 0
        y1=feval(func,t,p,c);
        func_calls = func_calls + 1;

        if (dp(j) < 0)		% backwards difference
            dydp(:,j) = (y1-y)./del(j);
        else			% central difference, additional func call
            p(j) = ps(j) - del(j);
	    dydp(:,j) = (y1-feval(func,t,p,c)) ./ (2 .* del(j));
            func_calls = func_calls + 1;
        end
     end

     p(j)=ps(j);		% restore p(j)

 end

% endfunction # ------------------------------------------------------ LM_DYDP


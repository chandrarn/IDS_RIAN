function [ sigma_p] = LM_Error(func,p,t,y_dat,dp,c)
%LM_Error Levenberg-Marquart final fit error function


 Npnt = length(y_dat);
 idx   = find(dp ~= 0);			% indices of the parameters to be fit
 Nfit = length(idx);

%delta_y = y_dat - feval(func,t,p_try,c);
delta_y = y_dat - feval(func,t,p,c) % using this version of the LM error
weight_sq = (Npnt-Nfit+1)/(delta_y'*delta_y) * ones(Npnt,1)
 
[alpha,~,~,~,~] = lm_matx(func,t,p,y_dat,weight_sq,dp,c)
covar = inv(alpha(idx,idx))
sigma_p = sqrt(diag(covar))

end


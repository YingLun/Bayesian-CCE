function [beta,res] = VAR_SUR(X,nlag,intercept)
%This function estimate a VAR(p) model by SUR-OLS.
%   INPUT:
%        X  = T*k matrix of observations
%     nlag  = p (scalar): number of lags
% intercept = logical, if intercept term is needed. (default = false)
%
%   OUTPUT:
%      beta = pN(+1)*k matrix of slope coefficients
%       res = (T-p)*k matrix of residuals


if nargin < 3
    intercept=false;
end
[no_obs, no_eqs] = size(X);

res = zeros(no_obs-nlag,no_eqs);
XLag = mlag(X,nlag);
if intercept
    XLag = [ones(no_obs,1),XLag];
    beta = zeros(nlag*no_eqs+1,no_eqs);
else
    beta = zeros(nlag*no_eqs,no_eqs);
end
XLag(1:nlag,:) = [];
for ii = 1:no_eqs
    y = X(nlag+1:no_obs,ii);
    beta(:,ii) = (XLag'*XLag)\(XLag'*y);
    res(:,ii) = y-XLag*beta(:,ii);
end

end


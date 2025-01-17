function [Gamma,R,Phi,Q,f] = StartVal(X,nlag,m,method)
%This function return the starting value of parameters.
%
%   INPUT:
%         X = T*k matrix of observations
%      nlag = p (scalar): number of lags
%         m = (scalar): number of factors
%    method = 1: OLS estimates
%           = 2: Null



if nargin == 3
    method = 1;
end

[T,k] = size(X);


switch method
    case 1
        % Eq.1
        [COEFF, SCORE, ~] = pca(X);
        f = SCORE(:,1:min(k,m));
        Gamma = COEFF(:,1:min(k,m));
        M = eye(T)-f*((f'*f)\f');
        e = M*X;
        R = e'*e/(T-k);
        % Eq.2
        [Phi,Res] = VAR_SUR(f,nlag);
        Phi = reshape(Phi,[m,m,nlag]);
        Q = cov(Res);
    case 2
        
    otherwise
        disp('Unknow method')
end


end


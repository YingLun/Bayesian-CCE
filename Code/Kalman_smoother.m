function [f_bac,Q_bac] = Kalman_smoother(f_tt,Q_tt,f_pred,Q_pred,Phi)
%This function carries out the Kalman smoother.
%(ONLY BACKWARD PASS, USE AFTER Kalman_filter.m)
%   INPUTS: 
%        f_tt = E[f(t)|x(1),...,x(t)]       (pm*T)
%        Q_tt = Cov(f_tt)                   (pm*pm*T)
%      f_pred = E[f(t)|x(1),...,x(t-1)]     (m*(T+1))
%      Q_pred = Cov(f_pred)                 (m*m*(T+1))
%         Phi = state coefficient matrix	(m*pm)
%
%   OUTPUTS:
%       f_bac = E[f(t)|x(1),...,x(T)]       (m*T)
%       Q_bac = Cov(f_bac)                  (m*m)

%% Housekeeping

[pm,T]          = size(f_tt);
[pm1,pm2,T1]    = size(Q_tt);
[m,Tp1]         = size(f_pred);
[m1,m2,Tp2]     = size(Q_pred);
[m3,pm3]        = size(Phi);

if (pm~=pm1)||(pm1~=pm2)||(pm2~=pm3)
    error('Kalman_smoother: Dimension pm not match.')
elseif (T~=T1)||(T1+1~=Tp1)||(Tp1~=Tp2)
    error('Kalman_smoother: Dimension T not match.')
elseif (m~=m1)||(m1~=m2)||(m2~=m3)
    error('Kalman_smoother: Dimension m not match.')
end

%% Initialization

f_bac	= zeros(pm,T);
Q_bac   = zeros(pm,pm,T);

f_bac(:,T)      = f_tt(:,T);
Q_bac(:,:,T)    = Q_tt(:,:,T);

disp('Kalman smoother: Started.')

for t=T-1:-1:1
    
    fprintf('Smoother %i of %i\n',t,T)
    Q_tt_tmp        = Q_tt(:,:,t);
    Q_pred_tmp      = Q_pred(:,:,t);
    Ct              = (Q_tt_tmp*Phi')/Q_pred_tmp;
    f_bac(:,t)      = f_tt(:,t)+Ct*(f_bac(1:m,t+1)-f_pred(:,t+1));
    Q_bac(:,:,t)    = Q_tt_tmp+Ct*(Q_bac(1:m,1:m,t+1)-Q_pred_tmp)*Ct';
end

f_bac = f_bac(1:m,:);
Q_bac = Q_bac(1:m,1:m,:);

disp('Kalman smoother: Finished.')

end


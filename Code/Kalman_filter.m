function [f_tt,Q_tt,f_pred,Q_pred] = Kalman_filter(X,Gamma,R,Phi,Q,ini_f,ini_Q)
%This function carries out the Kalman filter.
%   MODEL: (State-Space form)
%       observation: x(t) = Gamma*f(t)+e(t),     e(t)~N(0,R)
%             state: f(t) = Phi*f(t-1)+v(t),     v(t)~N(0,Q)
%
%   INPUT:
%           X = k*T observation matrix
%       Gamma = k*m coefficient matrix in observation equation
%           R = k*k covariance matrix of error terms in observation equation
%         Phi = m*m coefficient matrix in state equation
%           Q = m*m covariance matrix of state
%       ini_f = m*1 initial conditional state vector
%       ini_Q = m*m initial conditional state covariance
%
%   OUTPUT:
%        f_tt = E[f(t)|x(1),...,x(t)]
%        Q_tt = Cov(f_tt)
%      f_pred = E[f(t)|x(1),...,x(t-1)]
%      Q_pred = Cov(f_pred)

%% Housekeeping

[k,T] = size(X);
[kg,mg] = size(Gamma);
[kR1,kR2] = size(R);
[mP1,mP2] = size(Phi);
m = length(ini_f);
[mQ1,mQ2] = size(ini_Q);

if (k~=kg) || (kg~=kR1) || (kR1~=kR2)
    error('Kalman_filter: Dimension k not match')
elseif (mg~=mP1) || (mP1~=mP2) ||...
        (mP2~=m) || (m~=mQ1) || (mQ1~=mQ2)
    error('Kalman_filter: Dimension m not match')
end

disp('Dimension checked.')

%% Initialization

f_tt    = zeros(m,T);
Q_tt    = zeros(m,m,T);
f_pred  = zeros(m,T+1);
Q_pred  = zeros(m,m,T+1);

f_pred(:,1)     = ini_f;
Q_pred(:,:,1)   = ini_Q;

%% Filter

disp('Kalman filter: Started.')

for t = 1:T
    
    fprintf('Filter %i of %i\n',t,T)
    Q_tmp           = Q_pred(:,:,t);
    H_pred          = Gamma*Q_tmp*Gamma'+R;
    eta_pred        = X(:,t)-Gamma*f_pred(:,t);
    
    M_tmp           = Gamma'/H_pred;
    
    f_tt(:,t)       = f_pred(:,t)+Q_tmp*M_tmp*eta_pred;
    Q_tt(:,:,t)     = Q_tmp-Q_tmp*M_tmp*Gamma*Q_tmp;
    
    f_pred(:,t+1)   = Phi*f_tt(:,t);
    Q_pred(:,:,t+1) = Phi*Q_tt(:,:,t)*Phi'+Q;
    
end

disp('Kalman filter: Finished.')

end


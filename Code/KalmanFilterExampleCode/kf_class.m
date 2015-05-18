% function [output]=kf_class(input) 
function [output]=kf_class(y,A,B,sigma_e,sigma_u) 
% 4.2. Initialization
% -------------------------------------------------------------------------
n       = size(B,1);
T       = length(y);
Xtt1    = zeros(n,T);
Xtt     = zeros(n,T);
XtT     = zeros(n,T);

P_tt1   = zeros(n,n,T);
P_tt    = zeros(n,n,T);
P_tT    = zeros(n,n,T);

X_00    = zeros(n,1);
P_00    = 100*eye(n);

% 4.3. Kalman Filter
% -------------------------------------------------------------------------
for t=1:T
    % Prediction
    if t==1
        Xtt1(:,t)       = B*X_00;
        P_tt1(:,:,t)    = B*P_00*B' + sigma_u;
    else
        Xtt1(:,t)       = B*Xtt(:,t-1);
        P_tt1(:,:,t)    = B*P_tt(:,:,t-1)*B' + sigma_u;
    end
    % Updating
    Omega       = A*P_tt1(:,:,t)*A' + sigma_e;
    Kt          = P_tt1(:,:,t)*A'/Omega;
    ytilde      = y(1,t) - A*Xtt1(:,t);
    Xtt(:,t)    = Xtt1(:,t) + Kt*ytilde;
    P_tt(:,:,t) = P_tt1(:,:,t) - Kt*Omega*Kt';
    Kt_gain(t)  = Kt;
end

output.Kt_gain  = Kt_gain;
output.Xtt1     = Xtt1;
output.Xtt      = Xtt;
output.P_tt     = P_tt;
output.P_tt1    = P_tt1;

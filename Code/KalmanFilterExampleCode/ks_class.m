function [output]=ks_class(input,B) 
% function [output]=kf_class(input) 

T       = length(input.Xtt);
Xtt     = input.Xtt;
Xtt1    = input.Xtt1;
P_tt    = input.P_tt;
P_tt1   = input.P_tt1;

%% Kalman Smoother
% -------------------------------------------------------------------------
XtT(:,T)    = Xtt(:,T);
P_tT(:,:,T) = P_tt(:,:,T);

for t = T-1 : -1 : 1
    Jt          = P_tt(:,:,t)*B'/P_tt1(:,:,t+1);
    XtT(:,t)    = Xtt(:,t)    + Jt*(XtT(:,t+1)-Xtt1(:,t+1));
    P_tT(:,:,t) = P_tt(:,:,t) + Jt*(P_tT(:,:,t+1)-P_tt1(:,:,t+1))*Jt';
end

output.XtT  = XtT;
output.P_tT = P_tT;

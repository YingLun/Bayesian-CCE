%% Generate data
ARDL_Simulation;
X = [X1,X2];

%% Parameters
p = 3;          % No. of lags
m = 2;          % No. of factors

%% Initial values of (Theta,f)
[Gamma,R,Phi,Q_nu,f] = StartVal(X,p,m);

%% State-space form
Phi_tilde = [Phi';eye(m*(p-1)), zeros(m*(p-1),m)];
Q = [Q_nu,zeros(p,(m-1)*p);zeros((m-1)*p,m*p)];
f_tilde_full = vec(f);
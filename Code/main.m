% This is the main function.

%% Generate Data and Settings

ARDL_Simulation;
X = [X1,X2];

% Settings
Spat_Corr_X = false;    % FALSE if R is diagonal.
Rec_after   = 10;
Total_it    = 50;

% Parameters
p       = 3;          % No. of lags
m       = 1;          % No. of factors
[T,N]   = size(X1);
k       = 2;

%% Gibbs Sampling

% Initial values of (Theta,f)
[Gamma,R,Phi,Q_nu,f] = StartVal(X,p,m);

No_to_rec   = Total_it-Rec_after;
Gamma_rec   = zeros(N*k,m,No_to_rec);
R_rec       = zeros(N*k,N*k,No_to_rec);
fhat_rec    = zeros(T,m,No_to_rec);
Q_nv_rec    = zeros(m,m,No_to_rec);
Phi_rec     = zeros(m,m*p,No_to_rec);

for it=1:Total_it
    
    %% Step 1: Kalman Filter
    
    % Kalman Filter
    [Phi_til,Q] = companion(Phi,Q_nu);
    Gamma_til   = [Gamma,zeros(N*k,m*(p-1))];
    ini_f       = zeros(m*p,1);
    % ini_Q = zeros(m*p);
    % ini_Q(1:m,1:m) = eye(m);
    % ini_Q       = eye(m*p);
    ini_Q       = Q;
    
    [f_tt,Q_tt,f_til_pred,Q_til_pred] = Kalman_filter...
        (X',Gamma_til,R,Phi_til,Q,ini_f,ini_Q);
    
    f_pred  = f_til_pred(1:m,:);
    Q_pred  = Q_til_pred(1:m,1:m,:);
    
    [Ef,EQ] = Kalman_smoother(f_tt,Q_tt,f_pred,Q_pred,Phi_til(1:m,:));
    
    f_hat   = mvnrnd(Ef',EQ);
    
    %% Step 2: Densities of Theta(t+1)
    
    % Observation equation
    
    %OLS estaimtes and residuals
    Gamma_hat   = (f_hat'*f_hat)\(f_hat'*X);
    e_hat       = X-f*Gamma_hat;
    
    %Prior
    Gamma_0     = zeros(m,N*k);
    
    if Spat_Corr_X
        
        % R is not assumed to be diagonal.
        Omega_0     = eye(m);
        R_0         = cov(e_hat);
        
        R_bar       = R_0+e_hat'*e_hat+...
            Gamma_hat'/(Omega_0+inv(f_hat'*f_hat))*Gamma_hat;
        
        R           = iwishrnd(R_bar,N*k);
        
        Omega_bar   = inv(inv(Omega_0)+f_hat'*f_hat);
        Gamma_bar   = Omega_bar*(f_hat'*f_hat)*Gamma_hat;
        GammaVec    = mat2vec(Gamma_bar);
        Gamma       = vec2mat(mvnrnd(GammaVec',kron(R,Omega_bar))',m,N*k)';
        
    else
        % R is assumed to be diagonal.
        alpha       = 3;
        beta        = 0.001;
        M_0         = eye(m);
        
        alpha_star  = diag(alpha+e_hat'*e_hat+...
            Gamma_hat'/(inv(M_0)+inv(f_hat'*f_hat))*Gamma_hat);
        beta_star   = beta+T;
        
        M_i         = M_0+f_hat'*f_hat;
        R           = diag(1./gamrnd(alpha_star,1./beta_star));
        Gamma_bar   = M_i*(M_0\Gamma_0+f_hat'*f_hat*Gamma_hat);
        GammaVec    = mat2vec(Gamma_bar);
        Gamma       = vec2mat(mvnrnd(GammaVec,kron(R,inv(M_i)))',m,N*k)';
        
    end
    
    
    % State equation
    
    %OLS estaimtes and residuals
    f_lag           = mlag(f_hat,p);
    f_lag(1:p,:)    = [];
    f_hat(1:p,:)    = [];
    
    
    Phi_hat         = ((f_lag'*f_lag)\(f_lag'*f_hat))';
    v_hat           = f_hat-f_lag*Phi_hat';
    
    %Prior
    Phi_0           = zeros(m,m*p);
    Q_nu0           = cov(v_hat);
    Omega_f0        = inv(kron(diag(1:p),diag(diag(Q_nu0))));    
    
    %Posterior
    
    Q_bar       = Q_nu0+v_hat'*v_hat+...
        Phi_hat/(Omega_f0+inv(f_lag'*f_lag))*Phi_hat';
    
    Q_nu        = iwishrnd(Q_bar,m+T);
    
    Omega_bar   = inv(inv(Omega_f0)+f_lag'*f_lag);
    Phi_bar     = Omega_bar*(f_lag'*f_lag)*Phi_hat';
    PhiVec      = mat2vec(Phi_bar);
    Phi         = vec2mat(mvnrnd(PhiVec',kron(Q_nu,Omega_bar))',m,m*p);
    
    
    if it>Rec_after
        Gamma_rec(:,:,it-Rec_after)     = Gamma;
        R_rec(:,:,it-Rec_after)         = R;
        fhat_rec(:,:,it-Rec_after)      = f_hat;
        Q_nv_rec(:,:,it-Rec_after)      = Q_nu;
        Phi_rec(:,:,it-Rec_after)       = Phi;
    end
    
    Phi         = reshape(Phi,[m,m,p]);
    
end
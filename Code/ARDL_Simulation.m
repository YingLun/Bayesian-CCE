close all
clear, clc

%% Setting
N = 30;
T = 1000;
gamma_max = 0;
max_rho_e = 0;
X1 = zeros(T-500,N);
X2 = zeros(T-500,N);
Y = zeros(T-500,N);

%% Common Factor
rho_f = 0.6;
if rho_f==1
    sigma_zf = 0.1^2;
else
    sigma_zf = sqrt(1-rho_f^2);
end
zeta_f = randn(T,1)*sigma_zf;
f = filter(1, [1 -rho_f], zeta_f);
sigma_f = sigma_zf/sqrt(1-rho_f^2);


for ii = 1:N
    %%  x_1
    rho_v1 = rand(1)*0.95;
    sigma_zv1 = sqrt(1-rho_v1^2);
    zeta_v1 = randn(T,1)*sigma_zv1;
    v1 = filter(1, [1 -rho_v1], zeta_v1);
    Gamma_1 = randn(1)*0.2 + 1/sigma_f;
    x1 = Gamma_1*f+v1;

    %%  x_2
    rho_v2 = rand(1)*0.95;
    sigma_zv2 = sqrt(1-rho_v2^2);
    zeta_v2 = randn(T,1)*sigma_zv2;
    v2 = filter(1, [1 -rho_v2], zeta_v2);
    Gamma_2 = randn(1)*0.2 + 1/sigma_f;
    x2 = Gamma_2*f+v2;
    
    %% u
    rho_e = rand(1)*max_rho_e;
    sigma_e = sqrt(1-rho_e^2);
    zeta_e = randn(T,1)*sigma_e;
    epsilon = filter([1 -rho_e], 1, zeta_e);
    gamma_i = rand(1)*gamma_max;
    u = gamma_i*f + epsilon;
    
    %% y
    phi = 0.1 + (0.6-0.1).*rand(1);
    
    theta_1 = 1;
%     delta_1 = theta_1*(1-phi);
    w1 = rand(1);
    delta_10 = w1*theta_1*(1-phi);
    delta_11 = (1-w1)*theta_1*(1-phi);
    
    theta_2 = -1;
%     delta_2 = theta_2*(1-phi);
    w2 = rand(1);
    delta_20 = w2*theta_2*(1-phi);
    delta_21 = (1-w2)*theta_2*(1-phi);
    
    res_y = delta_10*x1(2:end)+delta_11*x1(1:end-1)+...
            delta_20*x2(2:end)+delta_21*x2(1:end-1)+u(2:end);
%     res_y = delta_1*x1+delta_2*x2+u;
    
    y = filter(1, [1 -phi], [0;res_y]);    
    
    X1(:,ii) = x1(501:end);
    X2(:,ii) = x2(501:end);
    Y(:,ii) = y(501:end);

end

clearvars -except X1 X2 Y
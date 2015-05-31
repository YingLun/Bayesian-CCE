function [Phi_til,Q] = companion(Phi,Q_nu)
%This function transform the VAR(p) model into its companion form.
%   MODEL: AR(p)-->AR(1)
%   f(t) = Phi(1)f(t-1)+..._Phi(p)f(t-p)+v(t),  v~N(0,Q_nu)
%   INPUT:
%       Phi = k*k*p array of k*k coefficients

%% Housekeeping
[m1,m2,p] = size(Phi);
if m1~=m2
    error('companion: Phi not square')
else
    m = m1;
end

[m1,m2] = size(Q_nu);
if (m1~=m)||(m2~=m)
    error('companion: Incorrect size of Q_nu')
end

%% Coefficient

Phi_til = zeros(m*p);
for ii=1:p
    Phi_til(1:m,(ii-1)*m+1:ii*m) = Phi(:,:,ii);
end
Phi_til(m+1:end,1:m*(p-1)) = eye(m*(p-1));


%% Covariance matrix
Q = zeros(m*p);
Q(1:m,1:m) = Q_nu;

end


function GammaVec = mat2vec(Gamma)
%This function vectorize Gamma.
%   INPUT
%       Gamma = m*N matrix
%   OUTPUT
%    GammaVec = mN*1 vector

[m,N] = size(Gamma);
GammaVec = zeros(N*m,1);
for ii=1:m
    GammaVec((ii-1)*N+1:ii*N) = Gamma(ii,:)';
end
end


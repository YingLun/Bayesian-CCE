function Gamma = vec2mat(GammaVec,m,N)
%This function vectorize Gamma.
%   INPUT
%    GammaVec = mN*1 vector   
%   OUTPUT
%       Gamma = m*N matrix

Gamma = zeros(m,N);
for ii=1:m
    Gamma(ii,:) = GammaVec((ii-1)*N+1:ii*N)';
end

end


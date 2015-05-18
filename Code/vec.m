function XVec = vec(X)
%This function vectorize X.

[T,k] = size(X);
XVec = zeros(T*k,1);
for ii=1:T
    XVec((ii-1)*k+1:ii*k) = X(ii,:)';
end
end


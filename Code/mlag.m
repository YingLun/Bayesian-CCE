function xlag = mlag(x,n)

if nargin ==1
    n = 1; % default value
end;

if nargin > 3
    error('mlag: Wrong # of input arguments');
end;

[nobs, nvar] = size(x);

xlag = zeros(nobs,nvar*n);
icnt = 0;
for ii=1:n
    for jj=1:nvar    
        xlag(ii+1:nobs,icnt+jj) = x(1:nobs-ii,jj);
    end;
    icnt = icnt+nvar;
end;

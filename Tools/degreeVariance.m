function [n,DV] = degreeVariance(gmt)

nmax = max(gmt(:,1));

A = zeros(nmax+1,nmax+1);
B = zeros(nmax+1,nmax+1);
k = 0;

for s=0:nmax
    
    i = s + 1;
    l = k + 1;
    k = l + nmax-s;
    
    A(i:end,i) = gmt(l:k,3);
    B(i:end,i) = gmt(l:k,4);
    
end

A = A.^2;
B = B.^2;

n = 0:nmax;
n = n';

DV = sum(A,2)+sum(B,2);
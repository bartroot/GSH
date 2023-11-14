function [Clm,Slm] = sc2zdenek(sc,lmax)
%
% Column vector of Clm and Slm coefficients to matrix of sc-format
% The maximum degree in sc can be larger than that in the vector.
%
%
% HOW
% [sc] = vec2sc(Clm,Slm,120)
%
% INPUT
% clm       vector of clm-coefficients
% slm       vector of slm-coefficients
%
% OUTPUT
% sc        coefficients in sc-format
%
% Wouter, January 27, 2007

lengthSHC = (lmax+1)/2 + (lmax+1)^2/2;

l1 = 0;
%sc = zeros(lmax+1,2*lmax+1);
Clm = zeros(lengthSHC,1);
Slm = zeros(lengthSHC,1);
temp = 1/(2*pi);
for ll = 0:lmax
    for mm = 0:ll
        l1 = l1+1;
        if l1<= length(Clm)            
            if (mm>0)
              Slm(l1,1) = -1*(-1)^mm.*sc(ll+1,lmax+1-mm)./sqrt(temp);
              Clm(l1,1) = (-1)^mm.*sc(ll+1,lmax+1+mm)./sqrt(temp);
            else
              Clm(l1,1) = (-1)^mm.*sc(ll+1,lmax+1+mm)./sqrt(0.5*temp);
            end
        else
        end
    end
end

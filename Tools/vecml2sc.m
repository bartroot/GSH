function [sc] = vecml2sc(Clm,Slm,lmax)
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

l1 = 0;
sc = zeros(lmax+1,2*lmax+1);
for mm = 0:lmax
    for ll = mm:lmax
        l1 = l1+1;
        if l1<= length(Clm)
            sc(ll+1,lmax+1+mm) = Clm(l1);
            if (mm>0)
              sc(ll+1,lmax+1-mm) = Slm(l1);
            end
        else
        end
    end
end

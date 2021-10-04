function [Clm,Slm,llvec,mmvec] = sc2vecml(sc,lmax)
%
% Column vector of Clm and Slm coefficients to matrix of sc-format
% The maximum degree in sc can be larger than that in the vector.
%
%
% HOW
% [Clm,Slm,llvec,mmvec] = sc2vecml(sc,120)
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

llvec = [];
mmvec = [];
Clm = zeros(   (lmax+2)*(lmax+1)/2  ,1  );
Slm = zeros(   (lmax+2)*(lmax+1)/2  ,1  );

for mm = 0:lmax
    for ll = mm:lmax
        l1 = l1+1;
        if l1<= length(Clm)
            Clm(l1) = sc(ll+1,lmax+1+mm);
            if (mm>0)
              Slm(l1) = sc(ll+1,lmax+1-mm);
            end
        else
        end
        llvec = [llvec ll];
        mmvec = [mmvec mm];
    end
end

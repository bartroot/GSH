function [V] = zdenek2V(VZ,maxdeg)
%
% function converts the SH coefficients from Zdenek's normalisation to GSH format

l = VZ(:,1);
m = VZ(:,2);

clm = VZ(:,3);
slm = VZ(:,4);

sc = zeros(maxdeg+1,2*maxdeg+1);
temp = 1/(2*pi);
for ll = 0:maxdeg
    for mm = 0:ll
        pp = find(l==ll&m==mm);        
        if ~isempty(pp)
            if pp<= length(clm)
                if (mm>0)
                  sc(ll+1,maxdeg+1+mm) =  (-1)^mm*sqrt(temp)*clm(pp);
                  sc(ll+1,maxdeg+1-mm) =  -1*(-1)^mm*sqrt(temp)*slm(pp);
                else
                  sc(ll+1,maxdeg+1+mm) =  (-1)^mm*sqrt(0.5*temp)*clm(pp);                  
                end
            end
        end        
    end  
end

[Clm,Slm,llvec,mmvec] = sc2vecml(sc,maxdeg);
V = [llvec' mmvec' Clm Slm];  
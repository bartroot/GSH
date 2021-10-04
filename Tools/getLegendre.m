function [setLeg] = getLegendre(SHbounds,Lat)
%
% This function computes the Legendre polynomial matrices used in the SH
% software package. Als the first and second derivative are computed for
% gradient computations
%
% input:  SHbounds = vector [nmin nmax], boundaries of the SH degree and order
%         Lat: latitude vector [degree] of the geogrid
%
% output: PMatrix: Legendre polynomials
%         DMatrix: first derivative
%         SMatrix: second derivative
%
% software routines used: visu2plm_ww
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nmin = SHbounds(1);
nmax = SHbounds(2);

CoLat = (90 - Lat(:,1))';
axl = (nmax)*(nmax+1)/2+nmax+1 - ((nmin-1)*(nmin)/2+nmin);

PnmF = zeros(axl,length(CoLat));
DnmF = zeros(axl,length(CoLat));
SnmF = zeros(axl,length(CoLat));

% Generate the Legendre function matrix for a particular theta
e = 1;
k = 1;
for s = 0:nmax

    b = e;
    
    % using wouters software to calculate polynomials
    [P,DP,SP] = Legendre_functions(nmin:nmax,s,CoLat);
    
    % Get the correct set (no zeros)
    newP = P';
    newDP = DP';
    newSP = SP';
    
    if s>nmin
        newP(1:k,:) = [];
        newDP(1:k,:) = [];
        newSP(1:k,:) = [];
        k = k + 1;
    end
    
    % file in temperal matrix
    PnmF(b:b + size(newP,1) - 1,:) =  newP;    
    DnmF(b:b + size(newDP,1) - 1,:) = newDP;
    SnmF(b:b + size(newSP,1) - 1,:) = newSP;
    
    e = b + size(newP,1);

end

% Making legendre setting
setLeg = struct();
setLeg.PMatrix = PnmF;
setLeg.DMatrix = DnmF;
setLeg.SMatrix = SnmF;

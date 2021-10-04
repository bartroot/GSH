function cs = sc2cs(field)

% SC2CS(FIELD) converts the rectangular (L+1)x(2L+1) matrix FIELD, containing
%       spherical harmonics coefficients in /S|C\ storage format into a 
%       square (L+1)x(L+1) matrix in |C\S| format.
%
% Nico Sneeuw
% Munich, 22/07/94

% uses none

[rows,cols] = size(field);
lmax = rows -1;
if cols ~= 2*lmax+1, error('Matrix dimensions must be (L+1)x(2L+1).'), end

c  = field(:,lmax+1:2*lmax+1);
s  = [zeros(lmax+1,1) field(:,1:lmax)];

cs = tril(c) + triu(rot90(s),1);

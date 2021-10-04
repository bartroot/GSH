function sc = cs2sc(field,backval)

% CS2SC(FIELD,backval) converts the square (L+1)x(L+1) matrix FIELD, containing
%       spherical harmonics coefficients in |C\S| storage format into a 
%       rectangular (L+1)x(2L+1) matrix in  /S|C\format.
%       The argument backval is optional and describes the matrix entries, 
%       where m > l. Default is 1e-20!

%----------------------------------------------------------------------------
% Nico Sneeuw, IAPG, TU-Munich                                  22/07/94
%----------------------------------------------------------------------------
% uses none
%
%----------------------------------------------------------------------------
% revision history
%   07/99, Matthias Weigelt - V5 adaptation, eliminiation of TRAPSTRIP
%----------------------------------------------------------------------------
% remarks
%
%----------------------------------------------------------------------------



if nargin == 1, backval = 1e-20; end
[rows,cols] = size(field);
lmax = rows -1;
if cols ~= rows, error('I expect a square matrix.'), end

c  = tril(field);
s  = rot90(triu(field,1),-1);

mask = backval*ones(lmax+1,2*lmax+1);
a = fliplr(triu(mask(:,1:cols-1)));
b = triu(mask(:,cols:2*lmax+1),1);
sc = [a b] + [s(:,2:lmax+1) c];

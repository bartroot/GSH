function [p, dp, ddp] = Legendre_functions(l,m,th)

% PLM Fully normalized associated Legendre functions for a selected order M
%
% HOW p       = plm(l,th)			- assumes M=0
%     p       = plm(l,m,th)
%     [p,dp]  = plm(l,m,th)
%     [p,ddp] = plm(l,m,th)
%
%
% IN  l  - degree (vector). Integer, but not necessarily monotonic.
%          For l < m a vector of zeros will be returned.
%     m  - order (scalar). If absent, m=0 is assumed.
%     th - co-latitude [deg] (vector)
% OUT p  - Matrix with Legendre functions. The matrix has length(TH) rows
%          and length(L) columns, unless L or TH is scalar. Then the output
%          vector follows the shape of respectively L or TH. 
%    dp  - Matrix with first derivative of Legendre functions. The matrix 
%          has length(TH) rows and length(L) columns, unless L or TH is 
%          scalar. Then the output vector follows the shape of respectively 
%          L or TH. 
%    ddp - Matrix with second (colatitude) derivative of Legendre functions. The matrix 
%          has length(TH) rows and length(L) columns, unless L or TH is 
%          scalar. Then the output vector follows the shape of respectively 
%          L or TH. 
% 
% Derivatives are based on GRafarend and Novak (2006) paper.
%-----------------------------------------------------------------------------
% Uses none
%-----------------------------------------------------------------------------
% Revision history:
%   - code based on visu2plm_ww from Nico Sneeuw and Wouter van der Wal
%   (23/11/2020) by Bart Root
%-----------------------------------------------------------------------------


% Input check for validatity of the program
if min(size(l)) ~= 1;  error('Degree l must be vector (or scalar)'); end
if any(rem(l,1) ~= 0); error('Vector l contains non-integers.'); end
if max(size(m)) ~= 1;  error('Order m must be scalar.'); end
if rem(m,1) ~=0;       error('Order m must be integer.'); end


% Preliminaries.
[lrow,lcol] = size(l);
[trow,tcol] = size(th);
lmax = max(l);
if lmax < m; error('Largest degree still smaller than order m.'); end
x    = cos(deg2rad(th(:)));
y    = sin(deg2rad(th(:)));
lvec = l(:)';					% l can be used now as running index.

% Recursive computation of the temporary matrix ptmp, containing the Legendre
% functions in its columns, with progressing degree l. The last column of
% ptmp will contain zeros, which is useful for assignments when l < m.
ptmp  = zeros(length(th),lmax-m+2);
dptmp = zeros(length(th),lmax-m+2); 
ddptmp = zeros(length(th),lmax-m+2);

%--------------------------------------------------------------------
% sectorial recursion: PM (non-recursive, though)
%--------------------------------------------------------------------
% WW: produces sqrt( (2n+1)/2n )
% Novak and Grafarend (2006), eq 64
% recursion is beta_n,n * beta_n-1,n-1 * ...
% sin(theta) * sin(theta) * ...
% * cos(theta) * P_0,0 (which is 1, dP_0,0 is 0)
% note that the term beta_n,n*cos(theta)*P_n-1,n-1 of Novak and
% Grafarend(2006) equation 72 is taken care of by the m.

% Calculate extra factor for derivatives of the Legendre functions.
if m == 0
   fac = 1;
else
   mm  = 2*(1:m);
   fac = sqrt(2*prod((mm+1)./mm));   % extra sqrt(2) because summation not over negative orders
end

% Calculate first column
ptmp(:,1) = fac*y.^m;                                      % The 1st column of ptmp.
dptmp(:,1) = m*fac*(y.^(m-1).*x);
ddptmp(:,1) = -m*fac*(y.^m) + m*(m-1)*fac*(y.^(m-2).*x.^2);


%--------------------------------------------------------------------
% l-recursion: P
%--------------------------------------------------------------------
for l = m+1:lmax
   col   = l - m + 1;			% points to the next column of ptmp
   root1 = sqrt( (2*l+1)*(2*l-1)/((l-m)*(l+m)) ) ;                      % beta_n,m (65) 
   root2 = sqrt( (2*l+1)*(l+m-1)*(l-m-1) / ( (2*l-3)*(l-m)*(l+m) ) );   % beta_n,m (65) * gamma_n,m (66)

   % recursion P
   if l == m+1
       ptmp(:,col) = root1 *x.*ptmp(:,col-1);
   else
       ptmp(:,col) = root1 *x.*ptmp(:,col-1) - root2 *ptmp(:,col-2);
   end
       
    % recursion DP
   if l == m+1
       dptmp(:,col) = root1 *(x.*dptmp(:,col-1)-y.*ptmp(:,col-1)); 
   else
       dptmp(:,col) = root1 *(x.*dptmp(:,col-1)-y.*ptmp(:,col-1)) - root2 *dptmp(:,col-2); 
   end

   % recursion DDP
   if l == m+1
       ddptmp(:,col) = root1 *(-x.*ptmp(:,col-1) -2*y.*dptmp(:,col-1) + x.*ddptmp(:,col-1) );
   else
       ddptmp(:,col) = root1 *(-x.*ptmp(:,col-1) -2*y.*dptmp(:,col-1) + x.*ddptmp(:,col-1) ) - root2*ddptmp(:,col-2); 
   end

end


% The Legendre functions have been computed. What remains to be done, is to
% extract the proper columns from ptmp, corresponding to the vector lvec. 
% If l or theta is scalar the output matrix p reduces to a vector. It should
% have the shape of respectively theta or l in that case.

lind       = find(lvec < m);			    % index into l < m
pcol       = lvec - m + 1;			        % index into columns of ptmp
pcol(lind) = (lmax-m+2)*ones(size(lind));	% Now l < m points to last col.

p          = ptmp(:,pcol);			        % proper column extraction 
dp         = dptmp(:,pcol);                 % proper column extraction
ddp        = ddptmp(:,pcol);                % proper column extraction

% Change order of vectors
if max(size(lvec))==1  && min(size(th))==1 && (trow == 1) 
    p = p'; 
    dp = dp';
    ddp = ddp';
end
if max(size(th))==1 && min(size(lvec))==1  && (lcol == 1) 
    p = p'; 
    dp = dp';
    ddp = ddp';
end

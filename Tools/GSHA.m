function cs = GSHA(f,lmax)

% GSHA global spherical harmonic analysis (wls)
%
% HOW cs = gsha(f,method,grid)		
%
% IN  f      - global field of size N*2N
%     lmax   - maximum degree of development
% OUT cs     - Clm & Slm in |C\S| format
%
%--------------------------------------------------------------------------
% uses Legendre_functions, SC2CS
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Grid definition
%--------------------------------------------------------------------------
[rows,cols] = size(f);

n     = rows;
dt    = 180 / n;
theta = (dt/2:dt:180)';
lam   = (dt/2:dt:360);           	% dt = dlam
   	
%--------------------------------------------------------------------------
% further diagnostics
%--------------------------------------------------------------------------
if lmax > n
    error('Maximum degree of development is higher than number of rows of input.')
end

if  lmax == n
    error('max. degree 180 is not possible for block grid, problem for m = 0')
end
%--------------------------------------------------------------------------
% Init.
%--------------------------------------------------------------------------
% L   = n;
L = lmax;
clm = zeros(L+1,L+1);
slm = zeros(L+1,L+1);

%--------------------------------------------------------------------------
% 1st step analysis: Am(theta) & Bm(theta)
%--------------------------------------------------------------------------
m   = 0:L;
c   = cos(lam'*m*pi/180);
s   = sin(lam'*m*pi/180);

% preserving the orthogonality 

c = c/rows; 
s = s/rows;

c(:,1)   = c(:,1)/2;			
s(:,L+1) = s(:,L+1)/2;          % not sure if this is needed? Need to take a look at the equations	
%c(:,L+1) = zeros(2*n,1); 		% this line is not needed otherwise Clamx,lamx is not computed.
s(:,1)   = zeros(2*n,1);   	
  
a = f * c;
b = f * s;

%--------------------------------------------------------------------------
% 2nd step analysis: Clm & Slm
%--------------------------------------------------------------------------
   
% predefine weights for the weighted least squares analysis
si = sin(theta*pi/180);
si = 2*si/sum(si);

% loop over the order of the Spherical Harmonics
for m = 0:L
  % construct the legendre polynomials
  p  = Legendre_functions(m:L,m,theta);
  % Select particular colum of the signal
  ai = a(:,m+1);
  bi = b(:,m+1);
  d   = 1:length(theta);
  pts = p' * sparse(d,d,si);
  
  % Estimate the coefficients
  clm(m+1:L+1,m+1) = (pts * p) \ pts * ai;
  slm(m+1:L+1,m+1) = (pts * p) \ pts * bi;  
end

%--------------------------------------------------------------------------
% Write the coefficients Clm & Slm in |C\S| format
%--------------------------------------------------------------------------
slm = fliplr(slm);
%sc  = [slm(:,1:L) clm];
cs  = sc2cs([slm(:,1:L) clm]);
cs  = cs(1:lmax+1,1:lmax+1);
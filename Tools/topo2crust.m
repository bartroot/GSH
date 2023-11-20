function [CM_interface,lon_CM,lat_CM] = topo2crust(h,N,CM_type,Model)

% this function initiates the global spherical harmonics analysis and synthesis 
% for general purpose | inserted functions are written by Bart Root, and 
% Weilun Qin from Delft University of Technology

%% define parameters and constants for Mars

Te = Model.Te;      			% Effective Elastic thickness of lithosphere[m]
T = Model.D_c;        			% Crustal thickness with zero topography [m] 
E = Model.E;       			% Young's modulus 
v = Model.v;       			% Poisson's ratio    
rho_m = Model.rho_m;   			% mantle density [kg/m^3]
rho_c = Model.rho_c;   			% crustal density [kg/m^3]
g = Model.GM./Model.Re^2;      		% gravity m/s^2
R = Model.Re;   			% radius in [m]

%% Calculate some variables
d_rho = rho_m - rho_c; 			% density difference [kg/m^3]f
D = E*Te^3/(12*(1-v^2)); 		% Flexural rigidity

disp(['The Flexural Rigidity used is: ' num2str(D./1e9,3) ' GPa m^3'])

% Airy local compensation
root = rho_c/d_rho*h;

%% Run GSHA (Analysis)
   
% Clm & Slm in |C\S| format [matrix]
cs = GSHA(root,N); 

% converts spherical harmonics coefficients in |C\S| storage format into 
% a rectangular (L+1)x(2L+1) matrix in  /S|C\format.
sc = cs2sc(cs);

% converts into SH degree, order and coefficients in matrix
[Clm,Slm,llvec,mmvec] = sc2vecml(sc,N);
V_root = [llvec' mmvec' Clm Slm];    

%% what type of crust will be initiated

if strcmp(CM_type,'Airy')

    V_CMroot = V_root.*1;

elseif strcmp(CM_type,'Infinite_Plate')

    n = V_root(:,1); 
    V_CMroot = V_root.*(1./(1+D/(d_rho)/g.*((2*n+1)./(2*R)).^4));

elseif strcmp(CM_type,'Thin_Shell')
    
    n = V_root(:,1); 
    AA = D/((d_rho)*g);
    BB = (n.^3.*(n+1).^3-4.*n.^2.*(n+1).^2+4.*n.*(n+1))./(n.*(n+1)-(1-v))./R^4;
    CC = 12*(1-v^2)/(R^2*Te^2);
    DD = (n.*(n+1)-2)./(n.*(n+1)-(1-v));
    V_CMroot = V_root.*(1./(1+AA.*(BB+CC.*DD)));

end

%% Run GSHS (Synthesis)

Clm_CM = V_CMroot(:,3); Slm_CM = V_CMroot(:,4);
sc_CM = vecml2sc(Clm_CM,Slm_CM,N);

% GSHS calculates a global spherical harmonic synthesis for any grid
% defined by lam and phi (each vectors). The radius must be scalar. The output
% is the disturbing potential and any derivative up to the fourth.

res = 180/(N+1);
lam = [res/2:res:360-res/2];         % lam   [n x 1]   longitude [deg]
th =  [res/2:res:180-res/2];         % th    [m x 1]   co-latitude [deg]      【90-lat】

% call the function
CM_root = GSHS(sc_CM,lam,th,N);
CM_interface = CM_root + T;

lon_CM = lam;
lat_CM = 90-th;


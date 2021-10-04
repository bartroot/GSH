function [data] = gravityModule(Lat,Lon,r,SHbounds,V,Re,GM)
%
% Calculating the the gravity field for the EIGENGL04C model
%
% input: Lat: latitude in degree [matrix]
%        Lon: longitude in degree [marix]
%        r: radial distance of computation surface (r=Re+h) [scalar]
%        SHbounds: [minSH max SH] example [0 5] or [2 7]
%        V : SH coefficients [same size as setleg files]
%        V format:  degree; order; Cnm; Snm
%        Re: radius of the Earth [meters]
%        GM: gravitational parameter of Earth 
%
% output: data: data structure of the gravity field
%
% programs used by the routine: 
%       - bsxfun
%       - getLegendre
%           - visu2plm_ww
%       
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the legendre polynomials

[setLeg] = getLegendre(SHbounds,Lat);

% Create the correct domain in the coefficients
nmin = SHbounds(1);
nmax = SHbounds(2);
V(V(:,1)>nmax|V(:,1)<nmin,:) = [];

% check is setLeg and V have the same length

if (size(setLeg.PMatrix,1)~=size(V,1))
    error('Length of Legendre polynomials should have the same length as the SH-coefficients')
end

% set data matrices

Pot   = zeros(size(Lon));
VecR  = zeros(size(Lon));
VecT  = zeros(size(Lon));
VecL  = zeros(size(Lon));
TenRR = zeros(size(Lon));
TenTT = zeros(size(Lon));
TenLL = zeros(size(Lon));
TenRT = zeros(size(Lon));
TenRL = zeros(size(Lon));
TenTL = zeros(size(Lon));

Trans1 = zeros(size(Lon));
Sint = zeros(size(Lon));

% Set some variables

fac  = GM/Re;
Sfac = GM/Re/Re;
Tfac = GM/Re/Re/Re;

lambda = deg2rad(Lon(1,:));
theta  = deg2rad(Lat(:,1))';
sint   = sin(pi/2-theta);
cost   = cos(pi/2-theta);

% start longitude for-loop, rangnge over all lon grid steps
lon = Lon(1,:);

for l = 1:length(lon)
    
    % Range related coefficients
    fac1 =  (Re./r).^(V(:,1)+1);
    Sfac1 = (Re./r).^(V(:,1)+2);
    Tfac1 = (Re./r).^(V(:,1)+3);

    %%%%%%%%%%%%%%%% Potential %%%%%%%%%%%%%%%%%%%%%%
    % Sumation of the Grace data potential

    pot = sum(bsxfun(@times,setLeg.PMatrix,((cos((V(:,2))*lambda(1,l)).*V(:,3)...
                                           + sin((V(:,2))*lambda(1,l)).*V(:,4))...
                                           .*fac1)),1);    

    %%%%%%%%%%%%%% Gradient vector %%%%%%%%%%%%%%%%%%
    % For the GRACE data files

    vecR =  sum(bsxfun(@times,setLeg.PMatrix,((cos((V(:,2))*lambda(1,l)).*V(:,3)...
                                             + sin((V(:,2))*lambda(1,l)).*V(:,4))...
                                             .*Sfac1.*(V(:,1)+1))),1);

    vecT =  sum(bsxfun(@times,setLeg.DMatrix,((cos((V(:,2))*lambda(1,l)).*V(:,3)...
                                             + sin((V(:,2))*lambda(1,l)).*V(:,4))...
                                             .*Sfac1)),1);

    vecL =  sum(bsxfun(@times,setLeg.PMatrix,((-(V(:,2)).*sin((V(:,2))*lambda(1,l)).*V(:,3)...
                                               +(V(:,2)).*cos((V(:,2))*lambda(1,l)).*V(:,4))...
                                               .*Sfac1)),1);

    %%%%%%%%%%%%%% Gradient tensor %%%%%%%%%%%%%%%%%%
    % Summation of the individual gradient elements

    tenrr =   sum(bsxfun(@times,setLeg.PMatrix,((cos((V(:,2))*lambda(1,l)).*V(:,3)...
                                               + sin((V(:,2))*lambda(1,l)).*V(:,4))...
                                            .*Tfac1.*(V(:,1)+1).*(V(:,1)+2))),1);

    tentt =   sum(bsxfun(@times,setLeg.SMatrix,((cos((V(:,2))*lambda(1,l)).*V(:,3)...
                                               + sin((V(:,2))*lambda(1,l)).*V(:,4))...
                                               .*Tfac1)),1);                

    tenll =   sum(bsxfun(@times,setLeg.PMatrix,((cos((V(:,2))*lambda(1,l)).*V(:,3)...
                                               + sin((V(:,2))*lambda(1,l)).*V(:,4))...
                                               .*Tfac1.*-(V(:,2)).^2)),1);

    tenrt =   sum(bsxfun(@times,setLeg.DMatrix,((cos((V(:,2))*lambda(1,l)).*V(:,3)...
                                               + sin((V(:,2))*lambda(1,l)).*V(:,4))...
                                               .*Tfac1.*(V(:,1)+1))),1);

    tenrl =   sum(bsxfun(@times,setLeg.PMatrix,((-(V(:,2)).*sin((V(:,2))*lambda(1,l)).*V(:,3)...
                                               +  (V(:,2)).*cos((V(:,2))*lambda(1,l)).*V(:,4))...
                                         .*Tfac1.*(V(:,1)+1))),1);

    tentl =   sum(bsxfun(@times,setLeg.DMatrix,((-(V(:,2)).*sin((V(:,2))*lambda(1,l)).*V(:,3)...
                                               +  (V(:,2)).*cos((V(:,2))*lambda(1,l)).*V(:,4))...
                                                 .*Tfac1)),1);   

    % Finalizing all the gravity data matrices

    Pot(:,l) = (fac*pot)';

    VecR(:,l) =  (Sfac*vecR)';
    VecT(:,l) = (-Sfac*vecT)';                    % added a minus sign
    VecL(:,l) =  (Sfac*vecL./sint)';

    TenRR(:,l) =  (Tfac*tenrr)';
    TenTT(:,l) =  (Tfac*tentt)';
    TenLL(:,l) =  (Tfac*tenll./sint./sint)';
    TenRT(:,l) =  (Tfac*tenrt)';                  % added a minus sign
    TenRL(:,l) = (-Tfac*tenrl./sint)';
    TenTL(:,l) = (-Tfac*tentl./sint)';            % added a minus sign
    
    % Transformation factors

    Trans1(:,l) = (cost./sint)';
    Sint(:,l) = sint';
    
end


% Construct data structure
data = struct();

% Gravity potential of topography
data.pot     = Pot;
data.grd.lon = Lon;
data.grd.lat = Lat;
data.grd.r   = r;

% Gravity vector of topography
data.vec.R = VecR;
data.vec.T = VecT;
data.vec.L = VecL;
data.vec.X = -VecT;
data.vec.Y = VecL;
data.vec.Z = -VecR;

% Gravity gradient derivatives - need to check these
data.ten.Trr = TenRR;
data.ten.Ttt = TenTT;
data.ten.Tll = TenLL;
data.ten.Trt = TenRT;
data.ten.Trl = TenRL;
data.ten.Ttl = TenTL;

% Gravity gradient Tensor in LNOU (see Zdenek's derivation)
data.ten.Txx = -VecR./r + TenTT;
data.ten.Tyy = -VecR./r + TenLL - Trans1.*VecT./r;
data.ten.Tzz = TenRR;
data.ten.Txy = (TenTL + Trans1.*VecL./r);
data.ten.Txz = -1.*(-VecT./r + TenRT);              % might need a multiplication with -1
data.ten.Tyz = -1.*(VecL./r - TenRL);               % might need a multiplication with -1
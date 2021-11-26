function [Gmap] = GreenMapping(depth,component,height,R,resolution,SHBounds)
%
% function computes Green's function mapping for fitting procedure
% relationships from MArtinec (2014)
%
% Input:   
%               depth    [meter]: depth of mass element wrt sea-level (R)
%               component [str} : Green's function type {Kr, KO, Krr, KrO, KOO}
%               height   [meter]: height of measurement wrt sea-level (R) [default = 0]
%               R        [meter]: Radius of sphere [default = 6378 km]
%               resolution [deg]: resolution of the map, equi-angular grid
%               SHBounds        : [minn maxn], sh resolution [default 0 220]
% 
% Output:
%               Gmap: normalised Green's function map 2D
%%%%%%

if nargin == 1
    error('Need at least a depth value for which the Green''s function is calculated and which Green''s function to compute')
elseif nargin < 3
    height = 0;
elseif nargin < 4
    height = 0;
    R = 6378000;    
elseif nargin < 5
    height = 0;
    R = 6378000;
    resolution = 1; % one degree
elseif nargin < 6
    height = 0;
    R = 6378000;
    resolution = 1; % one degree
    SHBounds = [0 220];
end

%% initialising 
j = (0:SHBounds(2));
m = 0;
size_map = 10; % in degrees (half size)

% angular distance between comp and integration point
phi = deg2rad((-size_map:resolution:size_map)); 
[X,Y] = meshgrid(phi,phi);

%Rm = ((sqrt(((cos(X)).^2+(cos(Y)).^2)./2)));

Rm = cos((sqrt((((X)).^2+((Y)).^2))));

r = R+height; 
rd = R-depth;
t = rd/r; % r/R

xR = (reshape(Rm,[],1));
%xP = (reshape(Rp,[],1));

%% get the Legendre polynomials, similar to Martinec (2014)

P = zeros(length(xR),SHBounds(2)+1);
DP = zeros(length(xR),SHBounds(2)+1);
DDP = zeros(length(xR),SHBounds(2)+1);
%dx = (x(2:end)+x(1:end-1))./2;
for n = j
   
   %P
   dum = legendre(n,(xR),'norm');
   pal0 = dum(1,:)';   
   P(:,n+1) = pal0;
   
   %dP/dx
%    pal0_mat = reshape(pal0,size(Rm,1),size(Rm,2));
%    [dx,dy] = gradient(pal0_mat);
%    dpal0_mat = sqrt(dx.^2+dy.^2);
%    dpal0 = reshape(dpal0_mat,[],1);
%    DP(:,n+1) = dpal0; 
   if n>0
        dpal0 = dum(2,:)';   
        DP(:,n+1) = dpal0;
   end
   
   %ddP/ddx
%    [ddx,ddy] = gradient(dpal0_mat);
%    ddpal0_mat = sqrt(ddx.^2+ddy.^2);
%    ddpal0 = reshape(ddpal0_mat,[],1);
%    DDP(:,n+1) = ddpal0;
    if n>1
        ddpal0 = dum(3,:)';   
        DDP(:,n+1) = ddpal0;
    end
    
end

%% removing long wavelengths TBD
% P(:,1:91) = 0;
% DP(:,1:91) = 0;
% DDP(:,1:91) = 0;
%

%% Green's functions: mapping of the function on a matrix 2mx2m, where m = size_map


%KrrP = sum((j+1).*(j+2).*t.^j.*P,2);
%KrOP = -sqrt(1-x'.^2).* sum((j+2).*t.^j.*DP,2);
%KOOP = 0.5.* (1-x'.^2) .* sum(t.^j.*DDP,2);

switch component
    case 'Kr'
        Gmap = - sum((j+1).*t.^j.*P,2);
    case 'KO'
        Gmap = - sum(sqrt(j.*(j+1)).*t.^j.*DP,2);
    case 'Krr'
        %Gmap = -1./sqrt(1+t.^2 - 2.*t.*(sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2))).^3 ...
        %                  + (3.*(1-t.*(sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2))).^2) ...
        %       ./(sqrt(1+t.^2 - 2.*t.*(sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2))).^5);
       %Gmap = sum((j+1).*(j+2).*t.^j.*P,2);
       
       Gmap = sum((j+1).*(j+2).*t.^j.*P.*sqrt((4.*pi)./(2.*j+1)),2);
    case 'KrO'
        %Gmap = -sqrt(1-sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2).^2).*(3*t.*(1-t.*sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2)))./(sqrt(1+t.^2 - 2.*t.*(sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2))).^5);        
        %Gmap = -sqrt(1-xR.^2) .* sum((j+2).*t.^j.*DP,2);
        
        %Gmap = sum((j+2).*sqrt(j.*(j+1)).*t.^j.*DP,2);
        Gmap = -sum((j+2).*sqrt(j.*(j+1)).*t.^j.*DP.*sqrt((4.*pi)./(2.*j+1)),2);
    case 'KOO'
        %Gmap = 0.5.*(1-sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2).^2).*(3*t.^2)./(sqrt(1+t.^2 - 2.*t.*sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2)).^5);
        %Gmap = 0.5.* (1-xR.^2) .* sum(t.^j.*DDP,2);
        
        %Gmap = 0.5.*sum(sqrt((j-1).*(j+1).*(j+2)).*t.^j.*DDP,2);
        Gmap = 0.5.*sum(sqrt((j-1).*(j+1).*(j+2)).*t.^j.*DDP.*sqrt((4.*pi)./(2.*j+1)),2);
    otherwise
        error('No such Green''s function component. Please use the following options: Kr, KO, Krr, KrO, or KOO.')
end

Gmap = reshape(Gmap,size(Rm,1),size(Rm,2));


        
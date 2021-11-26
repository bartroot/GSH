function [lampda,phi] = GreenFunctionv2(depth,component,height,R,SHBounds)
%
% function computes Green's function mapping for fitting procedure
% relationships from MArtinec (2014)
%
% Input:   
%               phi      [deg]  : lenght away from [point
%               depth    [meter]: depth of mass element wrt sea-level (R)
%               component [str} : Green's function type {Kr, KO, Krr, KrO, KOO}
%               height   [meter]: height of measurement wrt sea-level (R) [default = 0]
%               R        [meter]: Radius of sphere [default = 6378 km]
%               SHBounds        : [minn maxn], sh resolution [default 0 220]
% 
% Output:
%               Gvalue: normalised Green's function value
%%%%%%

if nargin < 2
    error('Need at least a depth value for which the Green''s function is calculated and which Green''s function to compute')
elseif nargin < 3
    height = 0;
    R = 6378000;
    SHBounds = [0 220];
elseif nargin < 4
    R = 6378000; 
    SHBounds = [0 220];
elseif nargin < 5       
    
    r = R+height;
    rd = R-depth;
    
    phi = (0:0.01:180);
    x = cos(deg2rad(phi));
    t = rd/r; % r/R

    g = sqrt(1+t.^2 - 2.*t.*x);

    %tensor mass-density Green's functions
    
    switch component
        case 'K'
            Gvalue = 1./(g);
        case 'Kr'
            Gvalue = -(1-t.*x)./(g.^3);
        case 'KO'
            Gvalue = sqrt(1-x.^2).*(t)./(g.^3);%- sum((j+1).*t.^j.*P,2);
        case 'Krr'           
           Gvalue = -1./g.^3 + (3.*(1-t.*x).^2)./(g.^5);         
        case 'KrO'
            Gvalue = -sqrt(1-x.^2).*(3*t.*(1-t.*x))./(g.^5);            
        case 'KOO'           
            Gvalue = 0.5.*(1-x.^2).*(3*t.^2)./(g.^5);
        otherwise
            error('No such Green''s function component. Please use the following options: Kr, KO, Krr, KrO, or KOO.')
    end    
else   

    %% initialising 
    j = (0:SHBounds(2));
    %m = 0;

    r = R+height; 
    rd = R-depth;
    t = rd/r; % r/R

    phi = 0:0.01:180;

    %% get the Legendre polynomials, similar to Martinec (2014)

    P = zeros(length(phi),SHBounds(2)+1);
    DP = zeros(length(phi),SHBounds(2)+1);
    DDP = zeros(length(phi),SHBounds(2)+1);
    %dx = (x(2:end)+x(1:end-1))./2;
    for n = j

       %P
       dum = legendre(n,cos(deg2rad(phi)),'norm');
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
            Gvalue = - sum((j+1).*t.^j.*P,2);
        case 'KO'
            Gvalue = - sum(sqrt(j.*(j+1)).*t.^j.*DP,2);
        case 'Krr'
            %Gmap = -1./sqrt(1+t.^2 - 2.*t.*(sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2))).^3 ...
            %                  + (3.*(1-t.*(sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2))).^2) ...
            %       ./(sqrt(1+t.^2 - 2.*t.*(sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2))).^5);
           %Gmap = sum((j+1).*(j+2).*t.^j.*P,2);

           Gvalue = sum((j+1).*(j+2).*t.^j.*P.*sqrt((4*pi)./(2.*j+1)),2);
           %Gvalue = sum((j+1).*(j+2).*t.^j.*P,2);
        case 'KrO'
            %Gmap = -sqrt(1-sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2).^2).*(3*t.*(1-t.*sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2)))./(sqrt(1+t.^2 - 2.*t.*(sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2))).^5);        
            %Gmap = -sqrt(1-xR.^2) .* sum((j+2).*t.^j.*DP,2);

            %Gmap = sum((j+2).*sqrt(j.*(j+1)).*t.^j.*DP,2);
            Gvalue = -sum((j+2).*sqrt(j.*(j+1)).*t.^j.*DP.*sqrt((4*pi)./(2.*j+1)),2);
            %Gvalue = sum((j+2).*sqrt(j.*(j+1)).*t.^j.*DP,2);
        case 'KOO'
            %Gmap = 0.5.*(1-sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2).^2).*(3*t.^2)./(sqrt(1+t.^2 - 2.*t.*sqrt((cos(deg2rad(X)).^2+cos(deg2rad(Y)).^2)./2)).^5);
            %Gmap = 0.5.* (1-xR.^2) .* sum(t.^j.*DDP,2);

            %Gmap = 0.5.*sum(sqrt((j-1).*(j+1).*(j+2)).*t.^j.*DDP,2);
            %Gvalue = 0.5.*sum(sqrt((j-1).*(j+1).*(j+2)).*t.^j.*DDP.*sqrt((4.*pi)./(2.*j+1)),2);
            Gvalue = 0.5.*sum(sqrt((j-1).*(j+1).*(j+2)).*t.^j.*DDP,2);
        otherwise
            error('No such Green''s function component. Please use the following options: Kr, KO, Krr, KrO, or KOO.')
    end

    %% correction value

    if ~strcmp(component,'KOO')
        KRR1 = -1/(sqrt(1+t^2-2*t).^3) + (3*(1-t)^2)./(sqrt(1+t^2-2*t)^5);
        Pdum = zeros(1,SHBounds(2)+1);
        for n = j
            dum = legendre(n,cos(deg2rad(0)),'norm');
            pal0 = dum(1,:)';   
            Pdum(1,n+1) = pal0;
        end
        Krr = sum((j+1).*(j+2).*t.^j.*P(1,:).*sqrt((4*pi)./(2.*j+1)),2);

        correctionfactor = Krr/KRR1;
    end
    if strcmp(component,'Krr')
        Gvalue = Gvalue/correctionfactor;
    end
    if strcmp(component,'KrO')
        Gvalue = Gvalue/correctionfactor;
    end

    %Total_G = abs(sum(Gvalue));
end

lampda = Gvalue;%./Total_G.*10; % times 10 because resolution is 0.1
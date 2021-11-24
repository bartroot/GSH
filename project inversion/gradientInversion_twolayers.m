function [density_inverted_crust,density_inverted_mantle,Hmatrix] = gradientInversion_v2(data,region,component,tikonov)
%
% This function will invert for density using the Green's functions as
% relation to the gravity data
%
% Input variables:
%                   data:           Similar to the data structure of the GSH package
%                   region:         Structure with information about the modelling region and bounds
%                   region.ax:      most northern latitude (deg)
%                   region.bx:      most southern latitude (deg)
%                   region.cx:      most western longitude (deg)
%                   region.dx:      most eastern longitude (deg)
%                   region.res:     resolution of inverted grid (deg)
%                   region.bound:   boundary around inverted region not used for solution (deg)
%                   components:     numerals for components 1:6 gradients, 7:9 vector, 10 geoid
%                   tikonov:        Tikonov regularisation parameter.
%
% Output variable:  
%                   density_inverted: Map of density specified by region.
%
% Functions used: GreenFunctionv2.m
%---------------------------------------------------------------------------

%component = 1;
%tikonov = 0; %[0 1e1 1e0 1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13 1e-14 1e-15 1e-16 1e-17 1e-18 1e-19 1e-20 1e-21 1e-22 1e-23 1e-24 1e-25 1e-26 -1e-27 1e-28 1e-29 1e-30 1e-31 1e-32 1e-33 1e-34 1e-35 1e-36 1e-37 1e-38 1e-39 1e-40 1e-41];

% Initializing parameters
height = data.height;
Re = data.Model.Re;
G = 6.6732e-11;
%SHbounds = data.SHbounds; % degrees of data

%% Modelling region

ax2 = region.ax;
bx2 = region.bx;

cx2 = region.cx;
dx2 = region.dx;

res = region.res;
bound = region.bound;

tbound = region.tbound;
mbound = region.mbound;
lbound = region.lbound;

depthc = -mean((tbound+mbound)./2,'all');
depthm = -mean((lbound+mbound)./2,'all');

% xxx1 = (ax2-bx2)/res + 1;
% xxx2 = (dx2-cx2)/res + 1;

%% Inversion region

ax = ax2+bound;
bx = bx2-bound;

cx = cx2 - bound;
dx = dx2 + bound;

xx1 = (ax-bx)/res + 1;
xx2 = (dx-cx)/res + 1;

%% Region of used observations

axo = ax;
bxo = bx;

cxo = cx;
dxo = dx;

% xxo1 = (axo-bxo)/res + 1;
% xxo2 = (dxo-cxo)/res + 1;

%% input model

% Downscale the WINTERC model to 2x2 degree
latLim = [bx ax res];
lonLim = [cx dx res];

lonD = lonLim(1):lonLim(3):lonLim(2);
latD = fliplr(latLim(1):latLim(3):latLim(2));
Lon = repmat(lonD,length(latD),1);
Lat = repmat(latD',1,length(lonD));

phi0 = deg2rad(90-Lat+res/2);
phi1 = deg2rad(90-Lat-res/2);

theta0 = deg2rad(Lon-res/2);
theta1 = deg2rad(Lon+res/2);

Vvolumec = (cos(phi0)-cos(phi1)).*(theta1-theta0).*(1./3.*(Re+tbound).^3-1./3.*(Re+mbound).^3);
Vvolumem = (cos(phi0)-cos(phi1)).*(theta1-theta0).*(1./3.*(Re+mbound).^3-1./3.*(Re+lbound).^3);

%%
new_Green = 1;
if new_Green == 1
    disp('Construct Green''s function')
    tic;
    [lamrr1,xrr1] =  GreenFunction(depthc,'Krr',height,Re);
    [lamrr2,xrr2] =  GreenFunction(depthm,'Krr',height,Re);
    toc
    [lamrO1,xrO1] =  GreenFunction(depthc,'KrO',height,Re);
    [lamrO2,xrO2] =  GreenFunction(depthm,'KrO',height,Re);
    toc
    [lamOO1,xOO1] =  GreenFunction(depthc,'KOO',height,Re);
    [lamOO2,xOO2] =  GreenFunction(depthm,'KOO',height,Re);
    toc
%     [lamr1,xr1] =  GreenFunction(depthc,'Kr',height,Re);
%     [lamr2,xr2] =  GreenFunction(depthm,'Kr',height,Re);
%     toc
%     [lamo1,xo1] =  GreenFunction(depthc,'KO',height,Re);
%     [lamo2,xo2] =  GreenFunction(depthm,'KO',height,Re);
%     toc
%     [lamK1,xK1] =  GreenFunction(depthc,'K',height,Re);
%     [lamK2,xK2] =  GreenFunction(depthm,'K',height,Re);
%     toc

    xrr1(isnan(lamrr1)) = [];
    lamrr1(isnan(lamrr1)) = [];
    xrO1(isnan(lamrO1)) = [];
    lamrO1(isnan(lamrO1)) = [];
    xOO1(isnan(lamOO1)) = [];
    lamOO1(isnan(lamOO1)) = [];
%     xr1(isnan(lamr1)) = [];
%     lamr1(isnan(lamr1)) = [];
%     xo1(isnan(lamo1)) = [];
%     lamo1(isnan(lamo1)) = [];
%     xK1(isnan(lamK1)) = [];
%     lamK1(isnan(lamK1)) = [];
    
    xrr2(isnan(lamrr2)) = [];
    lamrr2(isnan(lamrr2)) = [];
    xrO2(isnan(lamrO2)) = [];
    lamrO2(isnan(lamrO2)) = [];
    xOO2(isnan(lamOO2)) = [];
    lamOO2(isnan(lamOO2)) = [];
%     xr2(isnan(lamr2)) = [];
%     lamr2(isnan(lamr2)) = [];
%     xo2(isnan(lamo2)) = [];
%     lamo2(isnan(lamo2)) = [];
%     xK2(isnan(lamK2)) = [];
%     lamK2(isnan(lamK2)) = [];

    %save('GreensF_30_225_3120.mat','lamrr','lamrO','lamOO','xrr','xrO','xOO')
else
    %load('GreensF_30_225_3120.mat')
end          

%%
row_size = size(data.grd.lat,1)*size(data.grd.lon,2);
col_size = 2*size(Lat,1)*size(Lon,2);
hcol = size(Lat,1)*size(Lon,2);
H = zeros(row_size,col_size);
Hfull = zeros(col_size,col_size);
Hyfull = zeros(col_size,1);

%% deze moet aangepast worden!!!!!!!!!
LonO = (data.grd.lon);
LatO = (data.grd.lat); 

disp('Update the normal equations') 

for gg_num = component
    i = 0;
    for la = 1:size(LatO,1)
        if i>0;fprintf(repmat('\b',1,numel(S)+1));end
        S = sprintf(['Processing: ' num2str(floor(la/size(LatO,1)*100)) '% done!']); 
        fprintf(S);                                  

        COLAT1 = cos(deg2rad(90-LatO(la,1))).*cos(deg2rad(90-Lat)); 
        COLAT2 = sin(deg2rad(90-LatO(la,1))).*sin(deg2rad(90-Lat));

        for lo = 1:size(LonO,2)
            i = i + 1;

            %% Green's function map
            Lmap = COLAT1+COLAT2.*cos(abs(deg2rad(LonO(1,lo)-Lon)));                   
            if gg_num > 1 && gg_num < 6
                Hmap = 2*( mod(abs(deg2rad(LonO(lo)-Lon)), 2*pi) < pi ) - 1;
                cos_Theta = (cos(deg2rad(90-LatO(la)))-cos(deg2rad(90-Lat)).*Lmap)./(sin(deg2rad(90-Lat)).*sin(acos(Lmap)));
                cos_Theta(isnan(cos_Theta)) = 0;
                Amap = mod(real(Hmap.*acos(cos_Theta)),2*pi);

                if sum(sum(isnan(Amap))) > 0
                 error('Amap(isnan(Amap)) = 0;')
                end
            end           

           switch gg_num
               case 1
                    Gmap1 = interp1(cos(deg2rad(xrr1)),lamrr1,Lmap,'linear','extrap');
                    Gmap1(Lmap<cos(deg2rad(180))) = 0;
                    if sum(sum(isnan(Gmap1))) > 0
                        error('Nan values')
                    end
                    
                    Gmap2 = interp1(cos(deg2rad(xrr2)),lamrr2,Lmap,'linear','extrap');
                    Gmap2(Lmap<cos(deg2rad(180))) = 0;
                    if sum(sum(isnan(Gmap2))) > 0
                        error('Nan values')
                    end
                    
                    H(i,1:hcol) = reshape((G.*(1./(Re+height)^3).*Gmap1.*Vvolumec),1,[]);
                    H(i,hcol+1:end) = reshape((G.*(1./(Re+height)^3).*Gmap2.*Vvolumem),1,[]);                    
                    
                    if i==1                                 
                     y = reshape(data.ten.Tzz',[],1);                         
                    end                   
               case 2
                    Gmap1 = interp1(cos(deg2rad(xrO1)),lamrO1,Lmap,'linear','extrap');
                    Gmap1(Lmap<cos(deg2rad(180))) = 0;
                    if sum(sum(isnan(Gmap1))) > 0
                        error('Nan values')
                    end
                    Gmap1 = Gmap1.* 2.*cos(Amap);
                    
                    Gmap2 = interp1(cos(deg2rad(xrO2)),lamrO2,Lmap,'linear','extrap');
                    Gmap2(Lmap<cos(deg2rad(180))) = 0;
                    if sum(sum(isnan(Gmap2))) > 0
                        error('Nan values')
                    end
                    Gmap2 = Gmap2.* 2.*cos(Amap);
                    
                    H(i,1:hcol)          = reshape((G.*(1./(Re+height)^3).*Gmap1.*Vvolumec),1,[]);               
                    H(i,hcol+1:end)      = reshape((G.*(1./(Re+height)^3).*Gmap2.*Vvolumem),1,[]); 
                    
                    if i==1                                 
                     y = reshape(data.ten.Txz',[],1);
                    end                   
               case 3   
                    Gmap1 = interp1(cos(deg2rad(xrO1)),lamrO1,Lmap,'linear','extrap');
                    Gmap1(Lmap<cos(deg2rad(180))) = 0;
                    if sum(sum(isnan(Gmap1))) > 0
                        error('Nan values')
                    end
                    Gmap1 = Gmap1.*-2.*sin(Amap);
                    
                    Gmap2 = interp1(cos(deg2rad(xrO2)),lamrO2,Lmap,'linear','extrap');
                    Gmap2(Lmap<cos(deg2rad(180))) = 0;
                    if sum(sum(isnan(Gmap2))) > 0
                        error('Nan values')
                    end
                    Gmap2 = Gmap2.*-2.*sin(Amap);
                    
                    H(i,1:hcol)          = reshape((G.*(1./(Re+height)^3).*Gmap1.*Vvolumec),1,[]);
                    H(i,hcol+1:end)          = reshape((G.*(1./(Re+height)^3).*Gmap2.*Vvolumem),1,[]);
                    
                    if i==1                                 
                     y = reshape(data.ten.Tyz',[],1);                         
                    end                    
               case 4    
                    Gmap1 = interp1(cos(deg2rad(xOO1)),lamOO1,Lmap,'linear','extrap');
                    Gmap1(Lmap<cos(deg2rad(180))) = 0;
                    if sum(sum(isnan(Gmap1))) > 0
                        error('Nan values')
                    end
                    Gmap1 = Gmap1.*cos(2.*Amap);
                    
                    Gmap2 = interp1(cos(deg2rad(xOO2)),lamOO2,Lmap,'linear','extrap');
                    Gmap2(Lmap<cos(deg2rad(180))) = 0;
                    if sum(sum(isnan(Gmap2))) > 0
                        error('Nan values')
                    end
                    Gmap2 = Gmap2.*cos(2.*Amap);
                    
                    H(i,1:hcol)          = reshape((G.*(1./(Re+height)^3).*Gmap1.*Vvolumec),1,[]);
                    H(i,hcol+1:end)          = reshape((G.*(1./(Re+height)^3).*Gmap2.*Vvolumem),1,[]);
                    
                    if i==1                                 
                     y = reshape((data.ten.Txx-data.ten.Tyy)',[],1);                                                                 
                    end                    
               case 5
                    Gmap1 = interp1(cos(deg2rad(xOO1)),lamOO1,Lmap,'linear','extrap');
                    Gmap1(Lmap<cos(deg2rad(180))) = 0;
                    if sum(sum(isnan(Gmap1))) > 0
                        error('Nan values')
                    end
                    Gmap1 = Gmap1.*-2.*sin(2.*Amap);
                    
                    Gmap2 = interp1(cos(deg2rad(xOO2)),lamOO2,Lmap,'linear','extrap');
                    Gmap2(Lmap<cos(deg2rad(180))) = 0;
                    if sum(sum(isnan(Gmap2))) > 0
                        error('Nan values')
                    end
                    Gmap2 = Gmap2.*-2.*sin(2.*Amap);
                    
                    H(i,1:hcol)          = reshape((G.*(1./(Re+height)^3).*Gmap1.*Vvolumec),1,[]);
                    H(i,hcol+1:end)      = reshape((G.*(1./(Re+height)^3).*Gmap2.*Vvolumem),1,[]);
                    
                    if i==1                                 
                     y = reshape(data.ten.Txy',[],1);                         
                    end
               case 6
                    Gmap1 = interp1(cos(deg2rad(xrr1)),lamrr1,Lmap,'linear','extrap');
                    Gmap1(Lmap<cos(deg2rad(180))) = 0;
                    if sum(sum(isnan(Gmap1))) > 0
                        error('Nan values')
                    end
                    Gmap1 = Gmap1.*(-1./2);
                    
                    Gmap2 = interp1(cos(deg2rad(xrr2)),lamrr2,Lmap,'linear','extrap');
                    Gmap2(Lmap<cos(deg2rad(180))) = 0;
                    if sum(sum(isnan(Gmap2))) > 0
                        error('Nan values')
                    end
                    Gmap2 = Gmap2.*(-1./2);
                    
                    H(i,1:hcol)          = reshape((G.*(1./(Re+height)^3).*Gmap1.*Vvolumec),1,[]);
                    H(i,hcol+1:end)      = reshape((G.*(1./(Re+height)^3).*Gmap2.*Vvolumem),1,[]);
                    
                    if i==1                                 
                     y = reshape((data.ten.Txx+data.ten.Tyy)',[],1);                   
                    end                    
%                 case 7
%                     Gmap = interp1(cos(deg2rad(xr)),lamr,Lmap,'linear','extrap');
%                     Gmap(Lmap<cos(deg2rad(180))) = 0;
%                     if sum(sum(isnan(Gmap))) > 0
%                         error('Nan values')
%                     end                                 
%                     H(i,:)          = reshape((G.*(1./(Re+height)^2).*Gmap.*Vvolume),1,[]);
%                     if i==1                                 
%                      y = reshape((data.vec.Z)',[],1);                        
%                     end
%                 case 8
%                     Gmap = interp1(cos(deg2rad(xo)),lamo,Lmap,'linear','extrap');
%                     Gmap(Lmap<cos(deg2rad(180))) = 0;
%                     if sum(sum(isnan(Gmap))) > 0
%                         error('Nan values')
%                     end   
%                     Gmap = Gmap.*cos(Amap);
%                     H(i,:)          = reshape((G.*(1./(Re+height)^2).*Gmap.*Vvolume),1,[]);                    
%                     if i==1                                 
%                      y = reshape((data.vec.X)',[],1);                         
%                     end                    
%                 case 9
%                     Gmap = interp1(cos(deg2rad(xo)),lamo,Lmap,'linear','extrap');
%                     Gmap(Lmap<cos(deg2rad(180))) = 0;
%                     if sum(sum(isnan(Gmap))) > 0
%                         error('Nan values')
%                     end   
%                     Gmap = Gmap.*-1.*sin(Amap);
%                     H(i,:)          = reshape((G.*(1./(Re+height)^2).*Gmap.*Vvolume),1,[]);
%                     if i==1                                 
%                      y = reshape((data.vec.Y)',[],1);                        
%                     end                    
%                 case 10
%                     Gmap = interp1(cos(deg2rad(xK)),lamK,Lmap,'linear','extrap');
%                     Gmap(Lmap<cos(deg2rad(180))) = 0;
%                     if sum(sum(isnan(Gmap))) > 0
%                         error('Nan values')
%                     end                            
%                     H(i,:)          = reshape((G.*(1./(Re+height)).*Gmap.*Vvolume),1,[]);
%                     if i==1                                 
%                      y = reshape((data.pot)',[],1);                         
%                     end                    
               otherwise
                   error('no such setting')
           end


        end
        fprintf('\n')
    end  

%%     

    disp('Taking the transpose of the H matrix!')
    disp('...This can take a while, for resolutions smaller than 1 degree...')  
    
    HT = H'*H;
    Hy = H'*y; 
    
    disp(['Summating the updated normal equations for gradient component: ' num2str(gg_num)])          
    Hfull = Hfull + HT;
    Hyfull = Hyfull + Hy;

end

%% Estimating the density profile.
disp('Start the least-squares...')
disp('...This can take a while!')
deltaR = (Hfull+tikonov*eye(size(Hfull)))\Hyfull;

%P_dom = inv(Hfull+lambda*eye(size(Hfull)));
%% Results post-processing

density_inverted = reshape(deltaR,xx1,2*xx2);

density_inverted_crust = density_inverted(:,1:xx2);
density_inverted_mantle = density_inverted(:,xx2+1:end);

% For L_curve statistics
%AX_B = norm(Hfull*deltaR-Hyfull,2);
%B = norm(deltaR,2);
%AX_B2 = norm((Hfull+tikonov*eye(size(Hfull)))*deltaR-Hyfull,2);
%ax_b = norm(DM*deltaR-Y,2);



% Crustal and mantle inversion of WINTERC52 model
%
% We will go for the inverse method. Use the shape of the residual of the
% geoid potential between the model and GOCE data and multiply this with a
% factor and add a bias. Then recalculate the gravity and calculate a
% misfit. Plot this misfit for several values of the factor and bias.

clear;
close all;
clc;

warning off

HOME = pwd;

addpath([HOME '/Data']);
addpath([HOME '/Tools']);

%% settings

%create_new_slab = 1;                % creating new slab files (takes some time)


%directory = '/Users/bartroot/TUDelft/Research/3DEarth/Deep_Earth/Data/WINTERC/WINTERC54_MODEL/';
directory = '/home/bart/GSHA/Data/Synthetic/';

ExperimentalCaseNumber = 'E001';
model_name = '3layer_25km';
new_inital_model = 0;

height_fitting = 225000.0;
thick_lay = 25;

components = 1:6;

list_aaa = [0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1 1.5 2 2.5 3 3.5 4 4.5 5];
%list_aaa = [19 20 21 22 23 24 25];
llist = length(list_aaa); 

% inversion settings
Ltrunc = 2;    
max_ite = 1;

% for region
resolution = 0.5; 
region = struct();
region.res = resolution;
region.ax = 25+resolution/2;
region.bx = -19-resolution/2;
region.cx = 85-resolution/2;
region.dx = 135+resolution/2;
region.bound = 20;

% settings during the inversion
settings = struct();
settings.layers = [1 2 3];
    
%% %%%%% boundaries %%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Downloading input files and observations')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Downloading the boundaries')
[basement,Lon,Lat] = gmt2matrix(load([directory 'CCrust05.topo.gmt']));
%[mohoC] =            gmt2matrix(load('/Users/bartroot/PhD/Data/Leicester/High_resolution_CCrust/CCrust05.t15.gmt'));

WINTERC_MOHO_DATA = load([directory 'Global_Moho.lis']);
FWMOHO = scatteredInterpolant(WINTERC_MOHO_DATA(:,2),WINTERC_MOHO_DATA(:,3),-WINTERC_MOHO_DATA(:,4)./1e3,'linear','nearest');
moho = FWMOHO(Lon,Lat);

WINTERC_LAB_DATA = load([directory 'Global_LAB.lis']);
FWLAB = scatteredInterpolant(WINTERC_LAB_DATA(:,2),WINTERC_LAB_DATA(:,3),-WINTERC_LAB_DATA(:,4)./1e3,'linear','nearest');
lab = FWLAB(Lon,Lat);

clear WINTERC_LAB_DATA WINTERC_MOHO_DATA FWLAB FWMOHO

disp('boundaries downloaded')

%% download densities
disp('Downloading the densities')
load([directory 'Synthetic_densities.mat']);
disp('densities downloaded')

%% Download the observations

disp('Download observations coefficients')
load([directory 'Coefficients.Synthetic_model_' model_name '.mat']);
Vobs = VImodel;

% uncertainty of the observations
var_obs = 0.1E-9;

%% Compute the intial model

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Setup intial model')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%% Start complete process

disp('Start working on water layer...')                

Model = struct();

Model.number_of_layers = 1;
Model.name = 'WINTERC_waterlayer';

% Additional variables
Model.GM = 3.9860004415E14;
Model.Re_analyse = 6378136.30;
Model.Re = 6378136.30;
Model.geoid = 'none';
Model.nmax = 359;  
Model.correct_depth = 0;  

% % Top layer
Model.l1.bound = [directory 'CCrust05.top.gmt'];
Model.l1.dens  = [directory 'CCrust05.rhow.gmt'];
Model.l2.bound = [directory 'CCrust05.topo.gmt'];
% Global Spherical Harmonic Analysis
[Vwater] = model_SH_analysis(Model);    
disp('...water layer analysed')

Vinitial  = Vwater;

% rest of the initial model
disp('Setup density variation in initial model') 

if strcmp(ExperimentalCaseNumber,'E001')
    irho_crust = rho_crust - 0.03.*rho_crust;
    irho_litho = rho_litho - 0.03.*rho_litho;
    irho_mantl = rho_mantl - 0.03.*rho_mantl;
    
    varC = 2;
    varL = 2;
    varM = 2;
elseif strcmp(ExperimentalCaseNumber,'E002')
	RC = (rand(size(rho_crust))*2-1);
	RL = rand(size(rho_litho))*2-1;
	RM = rand(size(rho_mantl))*2-1;

	B = fspecial('gaussian');
	for rai = 1:5
		RC = filter2(B,RC);
		RL = filter2(B,RL);
		RM = filter2(B,RM);
	end
	varC = 100;
	varL = 250;
	varM = 100;

	irho_crust = rho_crust - RC*varC;
        irho_litho = rho_litho - RL*varL;
        irho_mantl = rho_mantl - RM*varM;

	save([directory 'Synthetic_densities_' model_name '_' ExperimentalCaseNumber '.mat'],'irho_crust','irho_litho','irho_mantl')
else
    error(['This Experimental case number does not exist: ' ExperimentalCaseNumber])
end

disp(['Experimental case number ' ExperimentalCaseNumber ' was set.'])

%%%% retrieve initial model coefficients
disp('Calculate gravity signal of inital model')

if new_inital_model == 1
    
    top_layer = -5:thick_lay:401;
    bot_layer = [top_layer(2:end) 401];

    for numl = 1:length(top_layer)
        disp(['Constructing initial layer number ' num2str(numl) '...'])

        ubound = -top_layer(numl);
        lbound = -bot_layer(numl);

        % LAB is shallower than 100 km
        upper_LAY = basement; 
        upper_LAY(basement>ubound) = ubound;                                                                                                                                       
        upper_LAY(basement<lbound) = lbound;

        mid_LAY = moho; 
        mid_LAY(moho>ubound) = ubound;                                                                                                                                       
        mid_LAY(moho<lbound) = lbound;

        lower_LAY = lab;
        lower_LAY(lab<lbound) = lbound;
        lower_LAY(lab>ubound) = ubound;

        % initiate Analysis file            
        Model.number_of_layers = 3;
        Model.name = 'Mantle layers';   

        Model.l1.bound = upper_LAY.*1e3;        
        Model.l1.dens  = irho_crust;

        Model.l2.bound = mid_LAY.*1e3; 
        Model.l2.dens = irho_litho ;%+ rho_iso;

        Model.l3.bound = lower_LAY.*1e3; 
        Model.l3.dens = irho_mantl;

        Model.l4.bound = lbound.*1e3;     
        % perform spherical harmonic analyses and synthesis       
        [Vlay] = model_SH_analysis(Model);

        % add to previous coefficients           
        Vinitial(:,3) = Vinitial(:,3) + Vlay(:,3);
        Vinitial(:,4) = Vinitial(:,4) + Vlay(:,4);    
    end    
    
    save([directory 'Coefficients_' model_name '_full_gradientfit_' ExperimentalCaseNumber '.mat'],'Vinitial')
else
    load([directory 'Coefficients_' model_name '_full_gradientfit_' ExperimentalCaseNumber '.mat'])
end
disp('Inital model is set')   

%% construct the a priory data

axa = region.ax+region.bound;
bxa = region.bx-region.bound;

cxa = region.cx - region.bound;
dxa = region.dx + region.bound;

xx1 = (axa-bxa)/region.res + 1;
xx2 = (dxa-cxa)/region.res + 1;
numdlay = xx1*xx2;
sizeAPR = numdlay*length(settings.layers);

settings.LAMBDA = zeros(sizeAPR,sizeAPR);
settings.dxapr = zeros(sizeAPR,1);

[acx,ccx,dum] = find(Lat==region.ax+region.bound&Lon==region.cx-region.bound);
[bcx,dcx,cum] = find(Lat==region.bx-region.bound&Lon==region.dx+region.bound);

for kk = 1:length(settings.layers)

   % Get data
    if settings.layers(kk) == 1
        dumLAM = diag((1./varC.^2).*ones(numdlay,1));
        dumapr = reshape(irho_crust(acx:bcx,ccx:dcx).*0,[],1);
    elseif settings.layers(kk) == 2
        dumLAM = diag((1./varL.^2).*ones(numdlay,1));
        dumapr = reshape(irho_litho(acx:bcx,ccx:dcx).*0,[],1);
    elseif settings.layers(kk) == 3
        dumLAM = diag((1./varM.^2).*ones(numdlay,1));
        dumapr = reshape(irho_mantl(acx:bcx,ccx:dcx).*0,[],1);
    else
        error('Unknown setting number in settings.layers')
    end
    
    settings.LAMBDA(1 + (numdlay*(kk-1)):numdlay*kk,1 + (numdlay*(kk-1)):numdlay*kk) = dumLAM;
    settings.dxapr(1 + (numdlay*(kk-1)):numdlay*kk,1) = dumapr;
                      
end

clear dumapr dumLAM
%% Weighting Matrix

settings.W = 1/(var_obs)^2.*eye(size(settings.LAMBDA));

%% Start the tikhonov iterations

VF = Vwater;

disp('Starting the inversion iterations')
for aaa = 1%:llist

    %tik_aaa = list_aaa(aaa);

    COMP = 'version3';%['synthv' num2str(tik_aaa) 'e20'];
    tikonov = NaN;%tik_aaa.*10.^(-(20));

    disp(COMP)

    %% boundaries of the regional inversion area

    disp('Set particular settings')

    region.tbound = basement(acx:bcx,ccx:dcx).*1e3;
    region.mbound = moho(acx:bcx,ccx:dcx).*1e3;
    region.lbound =  lab(acx:bcx,ccx:dcx).*1e3;    
    region.bbound =  -401.*ones(size(region.lbound)).*1e3;   
   
    %% Calculate the correction to the lithosphere such that the model fits the gravity data    

    %figure
    for ite = 1:max_ite       
        disp(['Gravity reduction iteration number: ' num2str(ite)])                
        if ite == 1          
                        
            % construct area for observed gravity field
            latLimo = [region.bx-region.bound region.ax+region.bound region.res];
            lonLimo = [region.cx-region.bound region.dx+region.bound region.res];

            SHboundso =  [Ltrunc 359];
            height =height_fitting;            

            % Observation coefficients
            %Vo = load('Data/XGM2016.txt');
            Vo = Vobs;
            Vo(Vo(:,1)>SHboundso(2),:) = [];                    
            Vo(3,3) = 0;

            % initialise densities
            drhoc = zeros(size(moho));
            drhol = zeros(size(moho));
            drhom = zeros(size(moho));
            
            % coefficients initial model
            Vi = Vinitial;            
            Vi(3,3) = 0;
            Vi(Vi(:,1)>SHboundso(2),:) = []; 

        end

        %% Calculating the difference
        disp(['Constructing the free-air anomaly of density model from iteration: ' num2str(ite-1)])

        Vres = Vo;
        Vres(:,3) = Vo(:,3) - Vi(:,3); 
        Vres(:,4) = Vo(:,4) - Vi(:,4);

        tic;
        [data_i] = model_SH_synthesis(lonLimo,latLimo,height,SHboundso,Vres,Model);               
        toc                                           

        %% Inversion of the density
        disp('Starting the inversion of the density with Green''s functaions') 
        [density_inverted_crust,density_inverted_litho,density_inverted_mantle,StatisticsV3] = gradientInversion_threelayers_V3(data_i,region,components,tikonov,settings);

        % inverted region
        [aix,cix,~] = find(data_i.grd.lat==region.ax&data_i.grd.lon==region.cx);
        [bix,dix,~] = find(data_i.grd.lat==region.bx&data_i.grd.lon==region.dx);                       

        %% Update densities: add densitities to initial model or previous iteration

        disp('Updating the drho matrix with the new reductionM values.') 
        reductionMc = zeros(size(drhoc));
        reductionMl = zeros(size(drhol));
        reductionMm = zeros(size(drhom));

        [amx,cmx,~] = find(Lat==region.ax&Lon==region.cx);
        [bmx,dmx,~] = find(Lat==region.bx&Lon==region.dx);
        reductionMc(amx:bmx,cmx:dmx) = (density_inverted_crust(bix:aix,cix:dix));
        reductionMl(amx:bmx,cmx:dmx) = (density_inverted_litho(bix:aix,cix:dix));
        reductionMm(amx:bmx,cmx:dmx) = (density_inverted_mantle(bix:aix,cix:dix));

        drhoc = drhoc + reductionMc;% - mean(mean(reductionM));        
        drhol = drhol + reductionMl;% - mean(mean(reductionM));
        drhom = drhom + reductionMm;% - mean(mean(reductionM));

	clear reductionMc redictionMl reductionMm
	% update the apriori knowledge
	%settings.LAMBDA = settings.LAMBDA + StatisticsV3.Pinv;
	%settings.dxapr = settings.dxapr + StatisticsV3.update_xapr;
	save([directory 'Stats_' num2str(ite) '.mat'],'StatisticsV3','settings')
        %% save the dhro variable to gmt file
        disp('Write drho file to memory')
        DRHOc = matrix2gmt(drhoc./1e3,Lon,Lat); 
        DRHOl = matrix2gmt(drhol./1e3,Lon,Lat); 
        DRHOm = matrix2gmt(drhom./1e3,Lon,Lat); 

        dlmwrite([directory 'WINTERC.drho_' num2str(ite) '_SH' num2str(Ltrunc) '_gradientfit_crust3_' ExperimentalCaseNumber '_' COMP '.gmt'],DRHOc,'delimiter','\t','precision','%15.8e')
        dlmwrite([directory 'WINTERC.drho_' num2str(ite) '_SH' num2str(Ltrunc) '_gradientfit_litho3_' ExperimentalCaseNumber '_' COMP '.gmt'],DRHOl,'delimiter','\t','precision','%15.8e')    
        dlmwrite([directory 'WINTERC.drho_' num2str(ite) '_SH' num2str(Ltrunc) '_gradientfit_mantle_' ExperimentalCaseNumber '_' COMP '.gmt'],DRHOm,'delimiter','\t','precision','%15.8e')    
        
	clear DRHOc DHROl DRHOm
        %% %%%%%%%%%%%%%% Update gravity signal of synthetic model %%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Calculating gravity signal of compensating masses')
                
        top_layer = -5:thick_lay:401;
        bot_layer = [top_layer(2:end) 401];

        for numl = 1:length(top_layer)
            disp(['Constructing layer number ' num2str(numl) '...'])

            ubound = -top_layer(numl);
            lbound = -bot_layer(numl);

            % LAB is shallower than 100 km
            upper_LAY = basement; 
            upper_LAY(basement>ubound) = ubound;                                                                                                                                       
            upper_LAY(basement<lbound) = lbound;

            mid_LAY = moho; 
            mid_LAY(moho>ubound) = ubound;                                                                                                                                       
            mid_LAY(moho<lbound) = lbound;

            lower_LAY = lab;
            lower_LAY(lab<lbound) = lbound;
            lower_LAY(lab>ubound) = ubound;

            % initiate Analysis file            
            Model.number_of_layers = 3;
            Model.name = 'Updates layers';   

            % crust
            Model.l1.bound = upper_LAY.*1e3;        
            Model.l1.dens  = irho_crust + drhoc;
            % lithosphere
            Model.l2.bound = mid_LAY.*1e3; 
            Model.l2.dens = irho_litho + drhol;
            % mantle
            Model.l3.bound = lower_LAY.*1e3; 
            Model.l3.dens = irho_mantl + drhom;
            Model.l4.bound = lbound.*1e3;     
            
            % perform spherical harmonic analyses and synthesis       
            [Vlay] = model_SH_analysis(Model);

            % add to previous coefficients           
            VF(:,3) = VF(:,3) + Vlay(:,3);
            VF(:,4) = VF(:,4) + Vlay(:,4);   

        end       

	clear Vlay upper_LAY mid_LAY lower_LAY
        disp('...Density model updated')
        
        Vi = VF;
        VF = Vwater;     
            
    end
    disp('Iterations finalised!')
    
    %% saving the coefficients
    disp('Saving Coefficients of latest model')
    V = Vi;
    
    save([ directory 'Coefficients_' model_name '_full_gradientfit_' ExperimentalCaseNumber '_' COMP '.mat'],'V')
    
    disp('End of the program!')

end

% main file for the complete GSH circle for a particular model
clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

% Model Construction

new_model = 1;

if new_model == 1

  % Construct new model
  
  Model = struct();
  
  Model.number_of_layers = 2;
  Model.name = 'two_layer_planet';
  
  % Additional variables
  Model.GM = 3.9860004415E14;
  Model.Re = 6378136.30;
  Model.geoid = 'none';
  Model.nmax = 179;     
  Model.correct_depth = 0;
  
  % Top layer
  Model.l1.bound = gmt2matrix(load([HOME '/Data/crust1.bd1.gmt'])).*1e3;  % meters with respect to reference sphere
  Model.l1.dens  = 2650;
  
  % Second layer
  Model.l2.bound = -50000;     % meters with respect to reference sphere
  Model.l2.dens  = 3400;	   % Density in kg/m3
  
  % Bottom bound
  Model.l3.bound = -100000;    % meters with respect to reference sphere
  
  % Save model in .mat file for use of the new software
  
  save([HOME '/Data/' Model.name '.mat'],'Model')

else
  % Load previous saved model

  model_name = 'two_layer_planet';
  load([HOME '/Data/' Model.name '.mat']);
end

%%%%%%%%%%%%%%%%%%% Computation area %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

latLim =    [-89.5 89.5 1];  % [deg] min latitude, max latitude, resolution latitude (preferable similar to latitude)
lonLim =    [-180 180 1];% [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
height =    10.0; % height of computation above spheroid
SHbounds =  [0 179]; % Truncation settings: lower limit, upper limit SH-coefficients used

%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

%% Global Spherical Harmonic Analysis 

tic;
[V_Model] = segment_2layer_model(Model.l1.bound,Model.l2.bound,Model.l3.bound,Model.l1.dens,Model.l2.dens,25000,Model)
toc

%% Global Spherical Harmonic Synthesis

tic;
[data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V,Model);
toc

%% Save data

DATE = datestr(now);
save(['Results/data_' Model.name '_' num2str(SHbounds(1)) '_' num2str(SHbounds(2)) '_' DATE '.mat'],'data','V','Model')

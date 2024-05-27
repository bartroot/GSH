function [data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V,Model)
% 
% This function is responsible for the SH synthesis of a given model.
%
% input:
%           - lonLim: longitude limits [min max resolution] in degree
%           - latLim: latitude limits  [min max resolution] in degree
%           - height: height of computation surface in meters [scalar/matrix]
%           - SHbounds: domain of order and degree SH coeff. [nmin nmax]
%           - V: full set of SHcoefficients constructed by model_SH_analysis.m
%           - Model: model structure constructed by inputModel.m contains Re and GM values
%
% output:   - data structure with gravity potential, gravity vector and gravity gradient tensor
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% construct the correct input variables

lon = lonLim(1):lonLim(3):lonLim(2);
lat = latLim(1):latLim(3):latLim(2);
    
Lon = repmat(lon,length(lat),1);
Lat = repmat(lat',1,length(lon));

r = Model.Re + height;

% make the Stokes coefficients independend of sorting

V = sortrows(V,2);

%% get potential field and others

if isscalar(r)
    [data] = gravityModule(Lat,Lon,r,SHbounds,V,Model.Re,Model.GM);
else
    [data] = gravityModule_full(Lat,Lon,r,SHbounds,V,Model.Re,Model.GM);
end

% saving the input values
data.latLim =    latLim; 
data.lonLim =    lonLim; 
data.height =    height;
data.SHbounds =  SHbounds;
data.Model = Model;

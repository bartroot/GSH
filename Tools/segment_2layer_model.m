function [V_Model] = segment_2layer_model(top_bound,middle_bound,bottom_bound,dens1,dens2,thickness_segment,Model)
%
% This function cuts a 2 layer model into segemnts to counteract the
% numerical instability of the GSH analysis code.
%
% input:            - top_bound             nxm matrix[in meters]
%                   - middle_bound          nxm matrix[in meters]
%                   - bottom_bound          scalar[in meters]
%                   - dens1                 nxm matrix, or scalar[in kg/m3]
%                   - dens2                 nxm matrix, or scalar[in kg/m3]
%                   - thickness_segment     scalar[in meters]
%                   - Model                 structure for Model.GM and Model.Re
%
% output:           - V_Model               nX4 vector containing SH coefficients
%
% Dependencies
%
% Needs to be used together with the GSH package
%
% written by Bart Root: 28 June 2023 for stdents PPI
% change record
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check input

% bottom bound needs to be a scalar
if ~isscalar(bottom_bound)
    error('bottom_bound needs to be a scalar!')
end
% thickness_segment needs to be a scalar
if ~isscalar(thickness_segment)
    error('thickness_segment needs to be a scalar!')
end
% does top segment encompass top segment
if thickness_segment<max(max(top_bound))
    top = 2*thickness_segment;
    if top<max(max(top_bound))
        top = 3*thickness_segment;
        if top<max(max(top_bound))
            top = 4*thickness_segment;
            if top<max(max(top_bound))
                error('top_bound is unable to be modeled.')
            end
        end
    end
else
    top = thickness_segment;
end

%% input parameters
bot = bottom_bound;
thick_lay = thickness_segment;
topo = top_bound;
moho = middle_bound;

%% check if first segment encompasses all topography

top_layer = top:-thick_lay:bot;
bot_layer = [top_layer(2:end) bot];
Model.number_of_layers = 2;

for numl = 1:length(top_layer)

    disp(['Constructing layer number ' num2str(numl) '...'])

    ubound = top_layer(numl);
    lbound = bot_layer(numl);

    % LAB is shallower than 100 km

    upper_LAY = topo; 
    upper_LAY(topo>ubound) = ubound;                                                                                                                                       
    upper_LAY(topo<lbound) = lbound;

    low_LAY = moho; 
    low_LAY(moho>ubound) = ubound;                                                                                                                                       
    low_LAY(moho<lbound) = lbound; 

    % layer1
    Model.l1.bound = upper_LAY;        
    Model.l1.dens  = dens1;

    % layer 2
    Model.l2.bound = low_LAY; 
    Model.l2.dens = dens2;

    % last bound
    Model.l3.bound = lbound; 

    % perform spherical harmonic analyses and synthesis       
    [Vlay] = model_SH_analysis(Model);

    % add to previous coefficients           
    if numl == 1
        V_Model = Vlay;
    else
        V_Model(:,3) = V_Model(:,3) + Vlay(:,3);
        V_Model(:,4) = V_Model(:,4) + Vlay(:,4);   
    end
end 
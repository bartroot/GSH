function [V] = model_SH_analysis(Model)
% 
% This function constructs the coefficients of the given model
%
% input: - Model structure (defined by inputModel.m)
%
% output: - V, SH coefficients of the complete model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Uses following functions:
%
%   - import_layer.m
%       - gmt2matrix.m
%   - layer_SH_analysis.m
%       - gsha_crust.m
%           - neumann.m (but not needed when wls is used)
%       - cs2sc.m
%       - sc2vecml.m
%       - geocradius.m (But this is an inert MATLAB function: aero toolbox)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Revision
%
% - by Bart Root, 30 june 2015: added fix for depth error, before going
%   into layer_SH_analysis.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the reference density
% 

%G = 6.67384e-11;      % Universal gravity constant
%G = 6.6732e-11;
%G = 6.67428e-11;
G = 6.673e-11;
Re3 = (Model.Re).^3;

rhoE = 3.*Model.GM./(4*pi*G*Re3);

% Default three binomial series terms (This can be modified after modifying layer_SH_analysis.m)

max_bin = 3;

% starting the loop over the different layers
for nlayer = 1:Model.number_of_layers
    
    layer_name = ['l' num2str(nlayer)];
    disp(['The layer number ' layer_name ' is starting'])
    
    % Costruct the coefficients for that particular layer
    [U,L,R] = import_layer(Model,nlayer);

    % costruct fix for depth error
    meanFix = max(max((U)));
    fixRe = Model.Re;
    Model.Re = Model.Re + meanFix;
    U = U - meanFix;
    L = L - meanFix;

    % start analysis fixed for depth error
    Vlayer = layer_SH_analysis(Model.nmax,Model.geoid,Model.Re,rhoE,max_bin,U,L,R);

    % transform back to reference sphere
    Vlayer(:,3) = (Model.Re./fixRe).^(Vlayer(:,1)+3) .* Vlayer(:,3);
    Vlayer(:,4) = (Model.Re./fixRe).^(Vlayer(:,1)+3) .* Vlayer(:,4);
    Model.Re = fixRe;      
    
    if isfield(Model.(layer_name), 'alpha')                
        [U,L,R,A] = import_layer(Model,nlayer);
        disp('Linear density gradient layer')
        
        [Vgradient] = gradient_SH_analysis(Model.nmax,Model.geoid,Model.Re,rhoE,max_bin,U,L,R,A);
        
        Vlayer(:,3) = Vlayer(:,3) + Vgradient(:,3);
        Vlayer(:,4) = Vlayer(:,4) + Vgradient(:,4);              
    end
        
    % Initialize at first iteration the full data set
    if nlayer==1
        V = zeros(size(Vlayer));
        V(:,1:2) = Vlayer(:,1:2);
    end
    
    % Summation of the different layer coefficients
    V(:,3:4) = V(:,3:4) + Vlayer(:,3:4);
    
end
    


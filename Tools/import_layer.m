function [U,L,R,A] = import_layer(Model,nlayer)
%
% This function imports the layer matrices

% name of the layer
layer_name = ['l' num2str(nlayer)];
layer_base = ['l' num2str(nlayer+1)];

% filenames
fupper = Model.(layer_name).bound;
flower = Model.(layer_base).bound;
fdens = Model.(layer_name).dens;

% check what class the input variables are
if (ischar(fupper) && ischar(flower))
    U = gmt2matrix(load(fupper)).*1e3;
    L = gmt2matrix(load(flower)).*1e3;
elseif (isnumeric(fupper) && isnumeric(flower))
    U = fupper;
    L = flower;
    if isscalar(U) && isscalar(L)
    else
        if isscalar(L)        
            L = flower.*ones(size(U));
        elseif isscalar(U)
            U = fupper.*ones(size(L));
        end
    end
elseif (ischar(fupper) && isnumeric(flower))
    U = gmt2matrix(load(fupper)).*1e3;
    L = flower;
    if isscalar(L)        
        L = flower.*ones(size(U));
    end
elseif (isnumeric(fupper) && ischar(flower))
    U = fupper;
    L = gmt2matrix(load(flower)).*1e3;
    if isscalar(U)        
        U = fupper.*ones(size(L));
    end
else
    error(['Could not resolve fupper and flower class & size:  ' class(fupper) ', ' class(flower)])
end

% Density file
if ischar(fdens)
    R = gmt2matrix(load(fdens)).*1e3;
elseif isnumeric(fdens)
    R = fdens;
    if isscalar(fdens)
        R = fdens.*ones(size(U));
    end
else
    error(['Could not resolve density class & size:  ' class(fdens)])
end

% both U and L scalar (only works if R isnot scalar)
if isscalar(R)
    error('All three input files are scalar, program needs info on size matrix!')
else
    if isscalar(U) && isscalar(L)
        U = fupper.*ones(size(R));
        L = flower.*ones(size(R));
    end
end   

% In case there is a linear density gradient in the model
if nargout == 4
    if isfield(Model.(layer_name), 'alpha')
        falpha = Model.(layer_name).alpha; 
        if ischar(falpha)
            A = gmt2matrix(load(falpha)).*1e3;
        elseif isnumeric(falpha)
            A = falpha;
            if isscalar(falpha)
                A = falpha.*ones(size(U));
            end
        else
            error(['Could not resolve linear density gradient class & size:  ' class(falpha)])
        end         
    else
        error('No specific file is listed for the linear gradient in the layer')
    end
end

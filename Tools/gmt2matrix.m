function [A,Lon,Lat] = gmt2matrix(c)
%
% This function converts the gmt imported files into matrices
%
% Written by Bart Root, 26-07-2012
%
% input:
%           - c(:,1) = longitude in degree
%           - c(:,2) = latitude in degree
%           - c(:,3) = value of the gmt file
%
%%%%%%%%%%%%%%% start of routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check input variable if matrix

%if ~ismatrix(c)
%    error('The input variabe must be a matrix');
%end

% Extract the latitude and longitude coordinates
res = abs(c(1,1)-c(2,1));
if res==0
    res = abs(c(1,2)-c(2,2));
end

lon = min(c(:,1)):res:max(c(:,1));
lat = min(c(:,2)):res:max(c(:,2));

A = zeros(length(lat),length(lon));
Lon = zeros(length(lat),length(lon));
Lat = zeros(length(lat),length(lon));

for n = 1:length(lat)
    for m =  1:length(lon)
        
        i = m + (n-1) * length(lon);
        
        A(n,m) = c(i,3);
        Lon(n,m) = c(i,1);
        Lat(n,m) = c(i,2);
        
    end
end

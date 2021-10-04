function [d] = matrix2gmt(A,Lon,Lat)
%
% This function converts the matrix into gmt based vectors
%
% Written by Bart Root, 26-07-2012
%
% input:
%           - A = value matrix
%           - Lon = longitude matrix in degree
%           - Lat = latitude matrix
%
%%%%%%%%%%%%%%% start of routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%

ll = length(Lat(:,1));
ww = length(Lat(1,:));
axl = size(A,1)*size(A,2);
d = zeros(axl,3);

for i = 1:ll
   
    bb = 1 + (i-1)*ww;
    ee = ww + (i-1)*ww;
    
    d(bb:ee,3) = A(i,:);
    d(bb:ee,1) = Lon(i,:); 
    d(bb:ee,2) = Lat(i,:);

end


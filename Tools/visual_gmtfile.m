function [A,Lon,Lat] = visual_gmtfile(filename,zunit,file_type)
%
% This function visualizes a .gmt file. 
%
% ONLY for inspection purposes

if nargin == 1
    file_type = 'block';
    zunit = '-';
end

if nargin == 2
    file_type = 'block';    
end

d = load(filename);

if strcmp(file_type,'block')
    [A,Lon,Lat] = gmt2matrix(d);
elseif strcmp(file_type,'gauss')
    [A,Lon,Lat] = gmt2matrix_gauss(d);
else
    error('File type must be string: block or gauss')
end

B = Europe_centered(A);

lon = Lon(1,:);
lats = Lat(:,1);
lon = lon - 180;

load coast;

figure;
hold on
imagesc(lon,lats,(B));c=colorbar; 
plot(long,lat,'k');
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
title(['Visualization of ' filename])
xlabel('LON degree')
ylabel('LAT degree')
ylabel(c,zunit)
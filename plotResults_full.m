% plot results
clear;
close all;
clc;

%% tutorial data
%%% insert output data file from Results here!!!%%%
load('Results/data_Crust10_crust_3_179_26-Mar-2019 13:14:45.mat')

%% plot different maps of the data
lon = data.grd.lon(1,:);
lats = data.grd.lat(:,1);

figure;
subplot(2,2,1)
imagesc(lon,lats,((data.pot)));c=colorbar; 
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Potential gravity field'])
ylabel(c,'m*m/s/s') 
set(gca,'YDir','normal')

subplot(2,2,2)
imagesc(lon,lats,((data.vec.Z)).*1e5);c=colorbar; 
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Z-component of gravity vector'])
ylabel(c,'mGal') 
set(gca,'YDir','normal')

subplot(2,2,3)
imagesc(lon,lats,((data.vec.X)).*1e5);c=colorbar; 
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['X-component of gravity vector (North-South)'])
ylabel(c,'mGal') 
set(gca,'YDir','normal')

subplot(2,2,4)
imagesc(lon,lats,((data.vec.Y)).*1e5);c=colorbar; 
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Y-component of gravity vector (East-West)'])
ylabel(c,'mGal') 
set(gca,'YDir','normal')

%% Tensor

figure;
subplot(3,3,1)
imagesc(lon,lats,((data.ten.Tzz).*1e9));c=colorbar;
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Tzz-component of gravity gradient tensor'])
ylabel(c,'Eotvos')
set(gca,'YDir','normal')

subplot(3,3,2)
imagesc(lon,lats,((data.ten.Txz).*1e9));c=colorbar;
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Txz-component of gravity gradient tensor'])
ylabel(c,'Eotvos') 
set(gca,'YDir','normal')

subplot(3,3,3)
imagesc(lon,lats,((data.ten.Tyz).*1e9));c=colorbar; 
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Tyz-component of gravity gradient tensor'])
ylabel(c,'Eotvos') 
set(gca,'YDir','normal')

subplot(3,3,5)
imagesc(lon,lats,((data.ten.Txx).*1e9));c=colorbar; 
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Txx-component of gravity gradient tensor'])
ylabel(c,'Eotvos') 
set(gca,'YDir','normal')

subplot(3,3,6)
imagesc(lon,lats,((data.ten.Txy).*1e9));c=colorbar;
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Txy-component of gravity gradient tensor'])
ylabel(c,'Eotvos') 
set(gca,'YDir','normal')

subplot(3,3,9)
imagesc(lon,lats,((data.ten.Tyy).*1e9));c=colorbar; 
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Tyy-component of gravity gradient tensor'])
ylabel(c,'Eotvos') 
set(gca,'YDir','normal')

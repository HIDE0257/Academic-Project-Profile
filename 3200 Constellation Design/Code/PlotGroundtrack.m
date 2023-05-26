%% Plots Groundtrack 
function PlotGroundtrack(lon,lat,F,V,targets,color,Nsc) 
v1 = V(F(targets,1),:); 
v2 = V(F(targets,2),:); 
v3 = V(F(targets,3),:); 
centroid = (v1+v2+v3)./3; 

[lon_target,lat_target,~] = cart2sph(centroid(:,1),centroid(:,2),centroid(:,3)); 
lon_target = rad2deg(lon_target); 
lon_target = 180 - lon_target;
lat_target = rad2deg(lat_target); 

cycle_start = [];
cycle_end = [];
n = 0;
for i = 1:length(lon)-1
    dlon = lon(i+1) - lon(i);
    if abs(dlon) >= 350
        n = n + 1;
        cycle_end(n,1) = i;
        cycle_start(n,1) = cycle_end(n) + 1;
    end
end
cycle_start = [1; cycle_start(1:end-1)];

figure(100)
for i = 1:length(cycle_start)
    p(1) = plot(lon(cycle_start(i):cycle_end(i)),lat(cycle_start(i):cycle_end(i)),'Color',color,'LineWidth',1, 'HandleVisibility','off'); hold on
    p(2) = plot(lon_target,lat_target,'*','Color','k', 'HandleVisibility','off'); hold on
    title(sprintf('Groundtacks of %d Satellites',Nsc))
    ylabel('Latitude [deg]')
    xlabel('Longtitude [deg]')
    axis([0 360 -90 90])
    grid on
end
end
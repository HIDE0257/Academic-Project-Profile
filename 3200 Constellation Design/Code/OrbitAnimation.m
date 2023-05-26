% ASEN 3200- Lab O1 Group Portion
% Date: 4/28/23
% Animation

clear;close all;clc;
load("New_Best_OE.mat")
v = VideoWriter('New_Best_OE_stationaryView.mp4','MPEG-4');

X_I = (zeros(length(t{1,1}),3)+1).*[-1 0 0];
  T_benu = 4.297461*3600;
  XI_body = ACI2body(X_I,T_benu,t{1,1});

open(v)

figure()
hold on
grid on
title("Constilation Propigation (ACI)")
patch('Faces',Facets,'Vertices',Verticies,'FaceColor','y')
v1 = Verticies(Facets(targets,1),:); 
v2 = Verticies(Facets(targets,2),:); 
v3 = Verticies(Facets(targets,3),:); 
centroid = (v1+v2+v3)./3; 
targetPlot = plot3(centroid(:,1),centroid(:,2),centroid(:,3),'.','Color','y','MarkerSize',10, 'HandleVisibility','off');
i = 1;
sc1 = plot3(X{1,1}(1:i,1),X{1,1}(1:i,2),X{1,1}(1:i,3),X{1,1}(i,1),X{1,1}(i,2),X{1,1}(i,3),"Color",'r');
sc2 = plot3(X{2,1}(1:i,1),X{2,1}(1:i,2),X{2,1}(1:i,3),X{2,1}(i,1),X{2,1}(i,2),X{2,1}(i,3),"Color",'b');
sc3 = plot3(X{3,1}(1:i,1),X{3,1}(1:i,2),X{3,1}(1:i,3),X{3,1}(i,1),X{3,1}(i,2),X{3,1}(i,3),"Color",'g');

for i=1:5:length(X{1,1})
    
plot3(X{1,1}(1:i,1),X{1,1}(1:i,2),X{1,1}(1:i,3),X{1,1}(i,1),X{1,1}(i,2),X{1,1}(i,3),"Color",'r');
plot3(X{2,1}(1:i,1),X{2,1}(1:i,2),X{2,1}(1:i,3),X{2,1}(i,1),X{2,1}(i,2),X{2,1}(i,3),"Color",'b');
plot3(X{3,1}(1:i,1),X{3,1}(1:i,2),X{3,1}(1:i,3),X{3,1}(i,1),X{3,1}(i,2),X{3,1}(i,3),"Color",'g');
h = gca;
plot3(h.XLim, [0,0], [0,0], 'k','LineWidth',2, 'HandleVisibility','off')
plot3([0,0], h.YLim, [0,0], 'k','LineWidth',2, 'HandleVisibility','off');
plot3([0,0], [0,0], h.ZLim, 'k','LineWidth',2, 'HandleVisibility','off');

%view(cross(X{2,1}(i,1:3),X{2,1}(i,4:6)))

view(XI_body(i,1:3))

%view(3)

axis equal

xlabel("X [km]",FontWeight="bold")

ylabel("Y [km]",FontWeight="bold")

zlabel("Z [km]",FontWeight="bold")
legend([targetPlot, sc1(1), sc2(1), sc3(1)],'Target Points', 'Spacecraft 1', 'Spacecraft 2', 'Spacecraft 3')

drawnow

frame = getframe(gcf);

writeVideo(v,frame);

end

hold off
close(v)

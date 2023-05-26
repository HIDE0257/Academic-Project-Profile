% Plot Orbits 
function plotOrbit(r_ACI,r_body,scn,Facets,Verticies,color,targets,viewMat,Nsc)
% This function takes orbitatl elements and produces an inital state vector
% to propigate the orbit 

r_ACI0 = r_ACI(1,:);
r_ACI_final = r_ACI(end,:);

v1 = Verticies(Facets(targets,1),:); 
v2 = Verticies(Facets(targets,2),:); 
v3 = Verticies(Facets(targets,3),:); 
centroid = (v1+v2+v3)./3; 

figure(1)
% figure
hold on
txt = ['Spacecraft ',num2str(scn)];
txt1 = ['Initial Position ', num2str(scn)];
txt2 = ['Final Position ', num2str(scn)];

title(Nsc+ " Constellation Orbits (ACI)")
plot3(r_ACI(:,1),r_ACI(:,2),r_ACI(:,3),"Color",color,'DisplayName',txt); hold on
plot3(centroid(:,1),centroid(:,2),centroid(:,3),'*','Color','g','MarkerSize',5, 'HandleVisibility','off');
plot3(r_ACI0(1),r_ACI0(2),r_ACI0(3),'.',"Color",'c','MarkerSize',25,'DisplayName',txt1);
plot3(r_ACI_final(1),r_ACI_final(2),r_ACI_final(3),'.','Color','k','MarkerSize',25,'DisplayName',txt2);
h = gca;
plot3(h.XLim, [0,0], [0,0], 'k','LineWidth',2, 'HandleVisibility','off')
plot3([0,0], h.YLim, [0,0], 'k','LineWidth',2, 'HandleVisibility','off');
plot3([0,0], [0,0], h.ZLim, 'k','LineWidth',2, 'HandleVisibility','off');
patch('Faces',Facets,'Vertices',Verticies,'FaceColor','y', 'HandleVisibility','off')
view(3)
xlabel("X_E [km]",FontWeight="bold")
ylabel("Y_E [km]",FontWeight="bold")
zlabel("Z_E [km]",FontWeight="bold")
grid on
hold off
legend show

figure(2)
% figure
hold on
title(Nsc+ " Constellation Orbits (Body)")
plot3(r_body(:,1),r_body(:,2),r_body(:,3),"Color",color,'DisplayName',txt');
h = gca;
plot3(h.XLim, [0,0], [0,0], 'k','LineWidth',2,'HandleVisibility','off')
plot3([0,0], h.YLim, [0,0], 'k','LineWidth',2, 'HandleVisibility','off');
plot3([0,0], [0,0], h.ZLim, 'k','LineWidth',2, 'HandleVisibility','off');
patch('Faces',Facets,'Vertices',Verticies,'FaceColor','y', 'HandleVisibility','off')
plot3(centroid(:,1),centroid(:,2),centroid(:,3),'*','Color','g','MarkerSize',5, 'HandleVisibility','off')
view(3)
xlabel("X_B [km]",FontWeight="bold")
ylabel("Y_B [km]",FontWeight="bold")
zlabel("Z_B [km]",FontWeight="bold")
legend show
grid on 
hold off

figure(3)
% figure
hold on
% title("Viewable Facet vs. Time (Spacecraft "+scn+")")
title("Viewable Facet vs. Time (" +Nsc+" Satellites)") 
scatter(viewMat(:,2),viewMat(:,1),[],rad2deg(viewMat(:,3)),'filled', 'HandleVisibility','off')
clrbar = colorbar;
clrbar.Label.String = 'Elevation Angle [deg]';
colormap jet
xlabel("Time [min]")
ylabel("Facet Number")
grid on
hold off


end
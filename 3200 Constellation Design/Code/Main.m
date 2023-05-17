
% Liam Mazzotta
% ASEN 3200- Lab O1 Group Portion
% Date: 4/28/23

clear;close all;clc;
tic
load("Target_list.mat",'targets') % 280 Points of intrest on Benu

%% Read in and Plot OBJ File
fileName ="Bennu.obj";
[Facets,Verticies] = read_obj(fileName,false,targets);

%% Main
%%% Sicence Value Objectives
% SV1: Number of spacecraft (minimize)
% SV2: Science Value, The cost function Jsv represents the amount of useful
% sciencetific information over the 1 week opservation peroid (maximize)
% Orbit constraints 0.300km <= rp <= ra <= 3km
% Make sure space craft do not colide with eachother or Benu

Nsc = 3; % Number of space craft
m = 1000/(1.05*Nsc); % mass [Kg]
A = (5-1/5*(Nsc-1))*1e-6; % Area [km^2]
mu_benu = 4.892e-9; %[km^3;/s^2]
t0 = 0; % Seconds
tf = 604800; % One Week in Seconds
T_benu = 4.297461*3600; % Peroid of Benu in seconds
XI = [-1 0 0]'; % Inertial Vector from stroid to the sun

%%% Convert Orbital Elements to initial state vector
% % mu, a, e, i, Omega, omega, theta
% X0_test = [0 -1 0 0 0 0.6994e-04];

%% 3 Spacecrafts 
% This is the best for 3 spadecrafts
a_sync = (T_benu*sqrt(mu_benu)/(2*pi))^(2/3);
e = [0 0.18 0.18];
a = [1.5*a_sync 1.5*a_sync 1.5*a_sync]; % [km]
incl = [pi/3 pi/6 -pi/4]; % rad] [70deg, 25deg -35deg]
Omega = [5*pi/12 pi/4 pi/4]; % [rad]
omega = [pi/2 pi/7 pi/7]; % [rad]
theta = [3*pi/2 pi 2*pi]; % [rad]

OE0 = [a; e; incl; Omega; omega; theta];
X0 = orbitalElements2X0(a,e,incl,Omega,omega,theta,mu_benu);

% %% 2 Spacecrafts
% a_sync = (T_benu*sqrt(mu_benu)/(2*pi))^(2/3);
% e = [0 0];
% a = [6*a_sync 6*a_sync]; % [km]
% incl = [0.6981 -0.6981]; % [rad] [70deg, 25deg -35deg]
% Omega = [pi/2 pi/2];    % [rad]
% omega = [pi/2 pi/2];    % [rad]
% theta = [0 pi];         % [rad]

%% 4 Spacecrafts
% a_sync = (T_benu*sqrt(mu_benu)/(2*pi))^(2/3);
% e = [0 0.18 0.18 0.18];
% a = [1.5*a_sync 1.5*a_sync 1.5*a_sync 1.5*a_sync]; % [km]
% incl = [pi/3 pi/6 -pi/4 pi/3]; % rad] [70deg, 25deg -35deg]
% Omega = [5*pi/12 pi/4 pi/4 -11*pi/36]; % [rad]
% omega = [pi/2 pi/7 pi/7 0]; % [rad]
% theta = [3*pi/2 pi 2*pi pi/2]; % [rad]
% X0 = orbitalElements2X0(a,e,incl,Omega,omega,theta,mu_benu);

% %% Best 
% % 4 Spacecrafts
% Nsc = 4;
% a_sync = (T_benu*sqrt(mu_benu)/(2*pi))^(2/3);
% e = [0 0 0.3 0.3];
% a = [1.5*a_sync 1.5*a_sync 1.8*a_sync 1.8*a_sync]; % [km]
% incl = [pi/3 -pi/3 2*pi/9 -2*pi/9]; % rad] 
% Omega = [pi/2 pi/2 pi/4 pi/4]; % [rad]
% omega = [pi/7 pi/7 pi/7 pi/7]; % [rad]
% theta = [pi/2 3*pi/2 pi 2*pi]; % [rad]
% X0 = orbitalElements2X0(a,e,incl,Omega,omega,theta,mu_benu);
%% 5 Spacecrafts
% a_sync = (T_benu*sqrt(mu_benu)/(2*pi))^(2/3);
% e = [0 0 0 0 0];
% a = [1.25*a_sync 1.25*a_sync 1.5*a_sync 1.5*a_sync 1.5*a_sync]; % [km]
% incl = [pi/3 -pi/3 pi/7 -pi/7 pi/2]; % rad] [70deg, 25deg -35deg]
% Omega = [pi/2 pi/2 pi/4 pi/4 0]; % [rad]
% omega = [pi/2 pi/2 pi/2 pi/2 0]; % [rad]
% theta = [pi/2 3*pi/2 pi 2*pi 0]; % [rad]
% X0 = orbitalElements2X0(a,e,incl,Omega,omega,theta,mu_benu);

color = ['r', 'b', 'g', 'm','y'];
for j = 1:Nsc
    [XOut,OEout,x,t_sc] = propagate_spacecraft(X0(j,:)',t0,tf,A,m);
    Final_OEout(j,1) = {OEout};
    t(j,1) = {t_sc};
    X(j,1) = {x};
    r_ACI{j,1} = X{j}(:,1:3);
    r_body(j,1) = {ACI2body(r_ACI{j},T_benu,t{j})};
    [observable{j,1}, elevationAngle{j,1}, cameraAngle{j,1},viewMat{j,1}] = check_view(t{j},r_body{j}, targets, Facets, Verticies,Nsc,j,T_benu);
    costMat = viewMat{j,1};
    GR = costMat(:,8);
    elevation = costMat(:,3);
    avg_elevation(j,1) = rad2deg(mean(elevation));
    H = costMat(:,6);
    JsvMat{j,1} = 0.007263./(GR.^2).*sin(elevation).*H;
    Jsv_vec(j,1) = sum([JsvMat{j,1}]);
    plotOrbit(r_ACI{j},r_body{j},j,Facets,Verticies,color(j),targets,viewMat{j},Nsc);
    % Groundtrack 
    for k = 1:length(t_sc)
        [lon(k,j),lat(k,j),dt_rev(k,j)] = groundtrack(t{j}(k),X{j}(k,:),T_benu);
    end
    PlotGroundtrack(lon(:,j),lat(:,j),Facets,Verticies,targets,color(j),Nsc);
end

AVG_elevation = mean(avg_elevation);

Jsv = sum(Jsv_vec);
fprintf("The science value Jsv: %E\n\n",Jsv)
fprintf('The average elevation angle: %.1f\n', AVG_elevation)
crashed_Sat = crashed_sat(Nsc,r_ACI);

%%% Functions
%% Read_OBJ
function [Facets,Vertices] =read_obj(fileName,plot,targets)
% Read in and Plot OBJ File

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name = "Benu"; color = 'y'; dataDelim2 = 'f';
% [facets,Verticies] =read_obj(fileName,Name,color,dataDelim2)
% INPUT
% fileName 
% Name (name of astroid)
% color (color of facets when ploted)
% dataDelim2 (vertical data differentiator)
%
% OUTPUT
%  3D plot of Astroid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Name = "Benu"; 
     color = 'y'; 
     dataDelim2 = 'f';

    FileID = fopen(fileName);
    Contents = textscan(FileID,'%c %f %f %f', 'Delimiter',{''},'CommentStyle',{'#'});
    Type = cell2mat(Contents(1,1));
    Index = find(Type== dataDelim2,1);
    Vector = [cell2mat(Contents(1,2)),cell2mat(Contents(1,3)),cell2mat(Contents(1,4))];
    Vertices = [Vector(1:Index-1,1) Vector(1:Index-1,2) Vector(1:Index-1,3)];
    Facets = [Vector(Index:end,1) Vector(Index:end,2) Vector(Index:end,3)];
    fclose(FileID);
    
    v1 = Vertices(Facets(targets,1),:); 
    v2 = Vertices(Facets(targets,2),:); 
    v3 = Vertices(Facets(targets,3),:); 
    centroid = (v1+v2+v3)./3; 
    if plot == true
        figure
        hold on
        title(Name)
        patch('Faces',Facets,'Vertices',Vertices,'FaceColor',color)
        plot(centroid,'*','r')
        xlabel("X [km]",FontWeight="bold")
        ylabel("Y [km]",FontWeight="bold")
        zlabel("Z [km]",FontWeight="bold")
        view(3)
        hold off
    end
end

%% OE to Initial State Vector
function X0 = orbitalElements2X0(a,e,i,Omega,omega,theta,mu)
n = 1;
for j = 1:length(a)
    % Deternime Position and Velocity in perifocial frame
    P = a(j)*(1-e(j)^2); 
    h = sqrt(mu*a(j)*(1-e(j)^2)); % Angular momentum
    
    r_PQW = [P*cos(theta(j))/(1+e(j)*cos(theta(j))); P.*sin(theta(j))/(1+e(j)*cos(theta(j))); 0];
    v_PQW = [-sqrt(mu/P)*sin(theta(j)); sqrt(mu/P)*(e(j)+cos(theta(j))); 0];
    
    % PQW to XYZ Rotation matrix
    T = [cos(Omega(j))*cos(omega(j))-sin(Omega(j))*sin(omega(j))*cos(i(j)), -cos(Omega(j))*sin(omega(j))-sin(Omega(j))*cos(i(j))*cos(omega(j)), sin(Omega(j))*sin(i(j)); sin(Omega(j))*cos(omega(j))+cos(Omega(j))*cos(i(j))*sin(omega(j)), -sin(Omega(j))*sin(omega(j))+cos(Omega(j))*cos(i(j))*cos(omega(j)), -cos(Omega(j))*sin(i(j)); sin(i(j))*sin(omega(j)), sin(i(j))*cos(omega(j)), cos(i(j))];    
    
    % Position and velocity in XYZ
    r_XYZ = T*r_PQW;
    v_PQW = T*v_PQW;
    
    % State Vector Output
    X0(n,1:6) = [r_XYZ;v_PQW]';
    n = n+1;
end
                
end

%% Spacecraft Propagate
function [XOut,OEout,X,t] = propagate_spacecraft(X0,t0,tf,A,m) 
% Function to propigate orbit with SRP pertibations, compute the orbit
% elements, position, and velocity in ACI at a given time t

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%  X0 - Spacecraft Cartesian in ACI frame at t0
%  t0 - scalar, Start of propigation (seconds)
%  tf - scalar, End of propigation (seconds)
%  A -  scalar, SRP area in km^2
%  m - scalar, spacecraft mass in kg
% OUTPUT:
% Xout = [x,y,z,x_dot,y_dot,z_dot]';
% OEout = [a,e,i,Omega,omega,theta]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu_benu = 4.892e-9; %[km^3;/s^2]

TSPAN = (t0:60:tf);  % one minute time step;
Opt = odeset('RelTol',10^-12,'AbsTol',10^-12);
[t,X] = ode45(@(t,X) twoBodyEOM(t,X,mu_benu,m,A),TSPAN,X0,Opt);

r_f = X(end,1:3)';
v_f = X(end,4:6)';
XOut = [r_f;v_f]; % Spacecraft cartesian state vector ACI frame at tf

epsilon = norm(v_f)^2/2-mu_benu/norm(r_f);
a = -mu_benu/(2*epsilon); % Semi-major axis of orbit

% Angular momentum
h = cross(r_f,v_f);

% Ecentricity of the orbit
e_vec = cross(v_f,h)/mu_benu-r_f/norm(r_f); 
e = norm(e_vec);

% Accending node Vector
n_vec = cross([0 0 1]',h);
n = norm(n_vec);

% Inclination (i)
%i = acosd(h(3)/norm(h)); % Deg
i = acos(h(3)/norm(h)); % Rad

% Right Ascension of Ascending Node (Ω)
if n_vec(2) < 0
   % Omega = 360 - acosd(n_vec(1)/n); % Deg
    Omega = 2*pi - acos(n_vec(1)/n); % Rad
else
    % Omega = acosd(n_vec(1)/n); % Deg
    Omega = acos(n_vec(1)/n); % Rad
end

% Argument of Perigee (ω)
if e_vec(3) < 0
    %omega = 360 - acosd(dot(n_vec,e_vec)/(n*e)); % Deg
    omega = 2*pi - acos(dot(n_vec,e_vec)/(n*e)); % Rad
else
    %omega = acosd(dot(n_vec,e_vec)/(n*e)); % Deg
    omega = acos(dot(n_vec,e_vec)/(n*e)); % Rad
end

if dot(r_f,v_f) < 0 
    %theta = 360 - acosd(dot(e_vec,r_f)/(e*norm(r_f)));% Deg
    theta = 2*pi - acos(dot(e_vec,r_f)/(e*norm(r_f)));% Rad
else
    %theta = acosd(dot(e_vec,r_f)/(e*norm(r_f))); % Deg
    theta = acos(dot(e_vec,r_f)/(e*norm(r_f))); % Rad
end

OEout = [a e i Omega omega theta]';

    function X_dot = twoBodyEOM(t,X,mu,m,A)

        r = X(1:3);
        r_dot = X(4:6);
        r_ddot = (-mu/(norm(r))^3).* r;
        
        % Acceleration due to solar radiation pressure
        P_sr = 4.57e-3;
        C_R = 1.2;
        a_SRP = -(P_sr*C_R*A/m)*[-1 0 0]';

        r_ddot = r_ddot + a_SRP; % Total acceleration

        X_dot = [r_dot; r_ddot];
    end

end

%% DCM_EtoB
function DCM_EtoB = T_EtoB(i,Omega,w)
% i = [rad] inclination
% Omega = [rad] right ascension
% w = [rad] argument of periapsis 

R1 = [cos(Omega) sin(Omega) 0; -sin(Omega) cos(Omega) 0; 0 0 1]; 
R2 = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)]; 
R3 = [cos(w) sin(w) 0; -sin(w) cos(w) 0; 0 0 1]; 

DCM_EtoB = [R3*R2*R1];

end

%% ACItoBody
function r_body = ACI2body(r_ACI,T_benu,t)
% This function rotates a vector in the ACI frame to the body frame

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% r_ACI - Position in ACI frame (Astroid Centerd Inertial Frame)
% T_benu - Peroid of rotation for Benu (seconds_)
% OUTPUT 
%  r_body- Position in body frame

Omega = 2*pi/T_benu; % rad/sec

    for i = 1:length(t)

        theta3 = Omega*t(i); % Rotation angle about benu z-axis

        % Rotation Matrix
        C = [   cos(theta3)     sin(theta3)     0;...
                -sin(theta3)    cos(theta3)     0;...
                     0              0           1       ];
         r_body(i,1:3) = r_ACI(i,1:3)*C;
         
    end
end

%% Check View
function [observable, elevationAngle, cameraAngle,viewMat] = check_view(t,r_body, targets, F, V, Nsc,scNumber,T_benu)
% Function to compute whether a space craft can observe a specified facet
% on the astroid surface

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% r - [x,y,z]' Spacecrafts cartesian position vector in body frame 3x1
% facetNumber - scalar, facet index
% F - matrix of verticies that form each face nx3
% V - matrix of vertex locations in implied body frame
% OUTPUT
% observable - int, 0 for unobservable, 1 for observable
% elevationAngle - scalar, elevation of spacecraft relative to facet plane in radians
% cameraAngle - scalar, angle of facet center relative to camera boresight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Observation requirments:
% local elevation angle must be at least 15 degrees
% center of the facet must be in the feild of view of the camera
% camera always points in the negative radial direction
    X_ACI = zeros(length(t),3) + [-1,0,0];
    X_I = ACI2body(X_ACI,T_benu,t);
    index = 1;
    for j = 1:length(t)
        r = [r_body(j,1) r_body(j,2) r_body(j,3)]';

         for k = 1:length(targets)
            facetNumber = targets(k);
            FacetIndex = F(facetNumber,:); % Verticies corosponding to facet index
            r_center = (sum(V(FacetIndex,:))/3)'; % Center pointiong vetor
            r_view = r-r_center; % Vector from camera to center of facet 
            d = norm(r_view)*1000;
            F_norm = cross(V(F(facetNumber,2),:)-V(F(facetNumber,1),:),V(F(facetNumber,3),:)-V(F(facetNumber,1),:));
            F_norm = F_norm/norm(F_norm);
            H = ceil(dot(X_I(j,:),F_norm));
            elevationAngle = asin(dot(r_view,F_norm)/norm(r_view));
            cameraAngle = acos(dot(r_view,r)/(norm(r_view)*norm(r)));
            FOV = pi/(9*Nsc);
            if elevationAngle >= pi/12 && cameraAngle <= FOV && H == 1
                observable = 1;
                GR = FOV/(2048-(1408/9)*(Nsc-1))*d;
                viewMat(index,1:8) = [facetNumber,t(j),elevationAngle,cameraAngle,FOV,H,d,GR]; 
                index = index+1;
            else
                observable = 0;
            end
           
        end
        
    end
end


%% Groundtrack 
function [lon,lat,dt_rev] = groundtrack(t,Xout,P) 
%     Input:  Xout = state vector at t [nx6] including 
%                     positions and velocties in the ACI frame.
%                     [km, km, km, km/min, km/min, km/min]
%             P = 1 revolution of Bennu spin [s] 
%             
%     Output: GT = groundtrack [Longtitude, Latitude] in [rad]

W_bennu = 2*pi/P;   % [rad/s] Angular velocity of Bennu

rev = floor(t/P);   % nth Revolution 
% Time for each revolution
if rev == 0
    dt_rev = t; 
elseif rev > 0 
    dt_rev = t - rev*P;
end

theta = W_bennu*dt_rev;     % [rad] Angle of Bennu relative to 
r_XYZ = Xout(1:3);          % [km] Satellite position in the ACI
R_XYZ = norm(r_XYZ);        % [km] Magnitude of satellite position in the ACI
r_xyz = T_EtoB(0,0,theta)*r_XYZ'; 
% R_xyz = norm(r_xyz); 

% lon = theta - acos(r_XYZ(1)/R_XYZ);     % [rad] Longtitude
% lat = atan(r_XYZ(3)*cos(lon)/r_XYZ(1)); % [rad] Latitude 
% lon = acos(r_xyz(1)/R_xyz) - theta;
% lat = atan(r_xyz(3)*cos(lon)/r_xyz(1));

[lon,lat,~] = cart2sph(r_xyz(1),r_xyz(2),r_xyz(3));
lon = rad2deg(lon);     % [deg] Longtitude 
lat = rad2deg(lat);     % [deg] Latitude [-90, 90]
if lon < 0 
    lon = 360 + lon;    % [deg] Longtitude [0, 360]
end

end

%% Crashed Satellite Check
function crashed = crashed_sat(Nsc,r_ACI)
% Inputs:   Nsc = Number of Satellites 
%           r_ACI = Satellite position in the ACI frame
% Outputs:  crashed = cell array {#ofCombination x 1}(#ofIndex x 3)[ith sat, jth sat, 1/0] 

% Ncomb = factorial(Nsc)/(factorial(2)*factorial(Nsc-2));
% crashed = zeros(length(r_ACI{1,1}),3);
comb = 2;
col_comb = 0;
for i = 1:Nsc-1
    for j = comb:Nsc
        col_comb = col_comb + 1;
        for k = 1:length(r_ACI{1,1})
            if r_ACI{i}(k,:) == r_ACI{j}(k,:)
                crashed{col_comb,1}(k,:) = {i, j, 1};
                fprintf('Our satellites (%d & %d) are crashed, baby!\n', i,j);
            else
                crashed{col_comb,1}(k,:) = {i, j, 0};
            end
        end
    end
    comb = comb + 1; 
end

end

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





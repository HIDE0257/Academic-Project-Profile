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



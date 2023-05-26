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

function [StreamFunction, Phi, Pressure, V, x, y] = Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N)
% [Input]: chord length, angle of attack, freestream velocity, freestream
% pressure, freestream density, number of panels

% [Output]: Streamlines, equipotential, pressure field, resultant velocity,
% x,y values

%% Define Domain
xmin = -0.5;
xmax = 6;
ymin = -1.5;
ymax = 1.5;

% Create meash over domatin using number of grid points specified
[x,y] = meshgrid(linspace(xmin,xmax,N),linspace(ymin,ymax,N));

r = @(x,y,x1) sqrt((x-x1).^2 + y.^2);           % Distance of each vortex
theta = @(x,y,x1) mod(atan2(y,(x - x1)),2*pi);  % Angle

dx = c/N;           %(m) Distance from a point to the other.
start_point = dx;   % Location of the 1st vortex
end_point = c - dx; % Location of the last vortex
point_x = linspace(start_point, end_point, N);   % Locations of each vortex

%% Streamline Function (Uniform Flow + Vortex)
% Initialize them to zero
Psi_V = 0;
Phi_V = 0;
u_vortex = 0;
v_vortex = 0;

% Uniform flow velocuty (u & v) in x and y
u = V_inf*cos(alpha);
v = V_inf*sin(alpha);
% Streamline & Potential for uniform flow
Psi_Uni = V_inf.*(y.*cos(alpha) - x.*sin(alpha));
Phi_Uni = V_inf.*(x.*cos(alpha) + y.*sin(alpha));

for i = 1:N
    %% Vortex Strength (Gamma(x))
    ratio_c = point_x(i)/c;     % (x/c) Ratio of each vortex location to the chord line
    gamma = 2*alpha*V_inf*sqrt((1 - ratio_c)/ratio_c);  % Vortex strength on the vortex sheet at given points
    Gamma = gamma*dx;           % Vortex strength
    
    %% StreamFunction for Vortex
    Psi_V = Psi_V + Gamma*log(r(x,y,point_x(i)))/(2*pi);
    
    %% Velocity Potential for Vortex
    Phi_V = Phi_V - Gamma*theta(x,y,point_x(i))/(2*pi);
    
    %% Pressure Distribution
    % Velocity in x and y for vortex
    u_vortex = u_vortex + Gamma.*sin(theta(x,y,point_x(i)))./(2.*pi.*r(x,y,point_x(i)));
    v_vortex = v_vortex - Gamma.*cos(theta(x,y,point_x(i)))./(2.*pi.*r(x,y,point_x(i)));
end

StreamFunction = Psi_Uni + Psi_V;   % Streamlines with all components combined
Phi = Phi_Uni + Phi_V;              % Velocity potential with all components combined
V_x = u + u_vortex;                 % (m/s) Total x-velocity
V_y = v + v_vortex;                 % (m/s) Total y-velocity
V = sqrt(V_x.^2 + V_y.^2);          % (m/s) Resultant velocity

q_inf = 0.5*rho_inf*V_inf^2;        % (Pa) Dinamic pressure
Pressure = p_inf + q_inf.*(1- (V./V_inf).^2);   % (Pa) Pressure distribution

end
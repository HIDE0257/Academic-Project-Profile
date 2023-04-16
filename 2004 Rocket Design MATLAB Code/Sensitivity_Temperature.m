clc; clear all
 
c = getConst(); % Call the getConst function to get constants

N = 50; 

l = 0.5;    % length of test stand
vx = 0; % initial x velocity of rocket
vy = 0; % initial y velocity of rocket
vz = 0; % initial z velocity of rocket
x0 = 0; % initial horizontal distance
y0 = 0; % initial horizontal distance
z0 = 0.25;  % initial vertical height

x = zeros(N,1);
y = zeros(N,1);
T_air = linspace(273.15, 311, N);

for i = 1:N
tspan = [0 20]; % time interval: 0 to 5 seconds
v_wind = normrnd(0,0)*0.45;
theta = normrnd(45,0); % initial angle of rocket
phi = normrnd(0,0);
%T_air = linspace(273.15, 311, N);
m_water0 = normrnd(1,0);
V_water0 = m_water0./c.rho_water;
V_air0 = c.Vbottle - V_water0;
P_atm = normrnd(83500,0);
CD = normrnd(0.2,0);
P0 = P_atm + c.P_gage;
m_air0 = P0.*V_air0./(c.R.*T_air(i));
m0 = c.m_bottle + m_air0 + m_water0;

val0 = [x0; y0; z0; vx; vy; vz; m0; V_air0; m_air0]; % initial variables for positions, velocities, mass, & volume of air
para = [v_wind, theta, phi, T_air(i), m_water0, P_atm, CD];
% ODE45 function [x; y; z; vx; vy; vz; m; V_air]
[t,Y] = ode45(@(t,Y) calcAccel(t,Y,l,para), tspan, val0);

x(i) = Y(end,1);
y(i) = Y(end,2);

figure(1)
hold on
grid on
plot(Y(:,1), Y(:,3));
title('Rocket Trajectory Varying Temperature of Water');
xlabel('Downrange Position (m)');
ylabel('Height (m)');

end 

figure(2)
grid on
plot(T_air,x,'bo');
title('Varying Temperature of Water');
xlabel('Temperature (K)');
ylabel('Downrange Distance (m)');
hold off 

% figure(2)
% hold on
% plot(x,y,'k.','markersize',6)
% axis equal; grid on; xlabel('x [m]'); ylabel('y [m]'); hold on;
%  
% % Calculate covariance matrix
% P = cov(x,y);
% mean_x = mean(x);
% mean_y = mean(y);
%  
% % Calculate the define the error ellipses
% n=100; % Number of points around ellipse
% p=0:pi/n:2*pi; % angles around a circle
%  
% [eigvec,eigval] = eig(P); % Compute eigen-stuff
% xy_vect = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
% x_vect = xy_vect(:,1);
% y_vect = xy_vect(:,2);
%  
% % Plot the error ellipses overlaid on the same figure
% plot(1*x_vect+mean_x, 1*y_vect+mean_y, 'b')
% plot(2*x_vect+mean_x, 2*y_vect+mean_y, 'g')
% plot(3*x_vect+mean_x, 3*y_vect+mean_y, 'r')
% title('Rocket Impact Locations & Error Ellipses');
% xlabel('Downrange Landing Distance (m)');
% ylabel('Crossrange Landing Distance (m)');
% hold off

 
% Sub function to calculate x,y,vx,vz,mass of the system, the air volume, mass of air
function A = calcAccel(t,y,l,para)
c = getConst(); % Call the getConst function to get all constants

mag_wind = para(1);
theta = para(2);
phi = para(3);
T_air = para(4);
m_water0 = para(5);
P_atm = para(6);
CD = para(7);
V_water0 = m_water0/c.rho_water;
V_air0 = c.Vbottle - V_water0;
P0 = P_atm + c.P_gage;
m_air0 = P0*V_air0/(c.R*T_air);


d = [y(1);y(2);y(3)];   % positions x and z
v = [y(4);y(5);y(6)];   % velocities x and z
v_wind = [mag_wind*cosd(phi); mag_wind*sind(phi); 0];
v = v - v_wind;         % relative velovity
m = y(7);               % mass of the system 
V_air = y(8);           % volume 
m_air = y(9);           % mass of air
 
M_v = sqrt(v(1)^2 + v(2)^2 + v(3)^2); % magnitude of velocity 
h = v./M_v; % the heading vector (unit vector for x and z)
 
% Set initial conditions for the heading vector. 
if sqrt(d(1)^2 + d(2)^2 + d(3)^2) <= l
    h = [cosd(theta); 0; sind(theta)];
end
 
At = pi*c.D_thrust^(2)/4; % Cross sectional area of the throat
Ab = pi*c.D_bottle^(2)/4; % Cross sectional area of the bottle 
 
P_end = P0*(V_air0/c.Vbottle)^c.gamma;  % Absolute pressure at the time all the water is expelled
P = P_end*(m_air/m_air0)^c.gamma;         % Absolute pressure at any time in Phase2
    
% Phase 1: Water expulsion
if  V_air < c.Vbottle % The mass of the system should be greater than m_phase1
    Pa = P0*(V_air0/V_air)^c.gamma; % absolute pressure at any time t
    delta_P = Pa - P_atm; 
    V_e = sqrt(2*delta_P/c.rho_water);  % pressure at the exit 
    
    mass = m; % the mass of the system
    fDrag = -0.5*c.rho_air*M_v^2*CD*Ab.*h; % drag forces (x and z) in Phase 1
    fThrust = 2*c.Cd*At*delta_P.*h;     % thrust forces (x ans z)in Phase 1
    fGrav = [0; 0; -mass*c.g];          % gravitational force in Phase 1
    v = [y(4);y(5);y(6)];   % store the values of velocity at times in Phase 1
    m = -c.Cd*At*sqrt(2*c.rho_water*delta_P); % the rate of change in mass of the system with time
    V_air = c.Cd*At*V_e;    % the rate of chang in voume of air with time
    m_air = 0;               % no change in mass of air in this phase
 
% Phase 2: Water Exhasted 
elseif P > P_atm
    P_cri = P*(2/(c.gamma+1))^(c.gamma/(c.gamma-1)); % Critical pressure
    rho_air_inside = m_air/c.Vbottle; % since air is compressible the air density varies at any time in Phase2 
    T = P/(rho_air_inside*c.R); % Accordingly, the temperature also varies at any time
    
    if P_cri > P_atm % if critical pressure is greater than the atm pressure
        % Assume the Mach number at the exit is 1 (M_e = 1: sonic flow)
        T_e = T*(2/(c.gamma+1)); % the exit temp 
        P_e = P_cri; % the exit pressure should be the same as the critical pressure
        rho_e = P_e/(c.R*T_e); % the exit air density
        V_e = sqrt(c.gamma*c.R*T_e); % the exit velocity 
    else  
        M_e = ((2/(c.gamma-1))*((P/P_atm)^((c.gamma-1)/c.gamma) - 1))^(0.5); % the Mach number at the exit
        T_e = T/(1+((c.gamma-1)*M_e^(2)/2)); % the exit temp
        rho_e =P_atm/(c.R*T_e); % the exit density
        P_e = P_atm; % the exit pressure should be equal to the atm pressure
        V_e = M_e*sqrt(c.gamma*c.R*T_e); % the exit velocity
    end
    
    mass = m; 
    fGrav = [0;0; -mass*c.g]; % gravitational force
    fDrag = -0.5*c.rho_air*M_v^2*CD*Ab.*h; % drage foces in x and z
    fThrust = (c.Cd*rho_e*At*V_e^2 + (P_atm - P_e)*At).*h; % thrust foces in x and z
    v = [y(4);y(5);y(6)]; % current velocities in x and z
    m = -c.Cd*rho_e*At*V_e; % the rate of change in mass of the system
    V_air = 0; % No change in volume of air inside the bottle
    m_air = -c.Cd*At*rho_e*V_e;
    
% Phase 3: Ballistic
else
    mass = m; % the mass of the system should be equal to mass of the bottle
    fGrav = [0;0;-mass*c.g]; % Gravitational force
    fDrag = -0.5*c.rho_air*M_v^2*CD*Ab.*h; % drag forces in x and z
    fThrust = [0;0;0]; % thrust forces in x and z
    v = [y(4);y(5);y(6)]; % current state of velocities in x and z
    m = 0; % no change in mass of the system 
    V_air = 0; % no change in volume of air
    m_air = 0; % no change in mass of air in this phase
    
    % When the rocket hits the ground, the rocket is at rest 
    if d(3) <= 0
        fGrav = [0;0;0];
        fDrag = [0;0;0];
        fThrust = [0;0;0];
        v = [0;0;0];
        m = 0.128;
        V_air = 0;
        m_air = 0;
    end
 
end
 
fnet = fDrag + fThrust + fGrav; % total force applied on the system
a = fnet./mass; % Acceleration
 
% state vector [vx; vz; ax; az; m; air volume; mass of air]
A = [v(1); v(2); v(3); a(1); a(2); a(3); m; V_air; m_air];
 
end
 
function c = getConst()
c.g = 9.81; % gravity (m/s^2)
c.Cd = 0.8; % discharge coefficient 
c.rho_air = 0.961; % ambient air density (kg/m^3)
c.Vbottle = 0.002; % volume of empty bottle (m^3)
c.gamma = 1.4; % ratio of the specific heat for air
c.rho_water = 1000;  % density of wter (kg/m^3)
c.D_thrust = 0.021; % diameter of throat (m)
c.D_bottle = 0.105; % diameter of bottle (m)
c.R = 287; % air gas constant (J/kgK)
c.m_bottle = 0.125; % mass of empty 2-litter bottle (kg)
c.P_gage = 40*6894.76; % initial gage pressure (Pa) (1 psi = 6894.76Pa) 
end
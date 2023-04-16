
clc; clear all
 
c = getConst(); % Call the getConst function to get constants

vx_wind = 0; 
vy_wind = 0;
vz_wind = 0;
v_wind = [vx_wind; vy_wind; vz_wind];

vx = 0; % initial x velocity of rocket
vy = 0; % initial y velocity of rocket
vz = 0; % initial z velocity of rocket
x0 = 0; % initial horizontal distance
y0 = 0; % initial horizontal distance
z0 = 0.25; % initial vertical height
theta = 47; % initial angle of rocket
l = 0.5; % length of test stand
 
tspan = [0 20]; % time interval: 0 to 5 seconds
val0 = [x0; y0; z0; vx; vy; vz; c.m0; c.V_air0; c.m_air0]; % initial variables for positions, velocities, mass, & volume of air
 
% ODE45 function [x; z; vx; vz; m; V_air]
[t,y] = ode45(@(t,y) calcAccel(t,y,theta,l,v_wind), tspan, val0);

z_max = max(y(:,3)); % Find the maximum height of the rocket
fprintf('Max height (m): %.2f\n', z_max);
x_max = max(y(:,1)); % Find the maximum displacement of the rocket
fprintf('Max distance (m): %.2f\n', x_max);
d_ang = atand(y(end,2)/y(end,1));
fprintf('Drifted Angle (deg): %.1f\n', d_ang);

% Plot x vs z
figure(1)
plot3(y(:,1),y(:,2),y(:,3),'-')
title('Trajectory of the Rocket', 'FontSize', 20);
xlabel('Displacement (m)','FontSize', 15);
ylabel('Drifted (m)', 'FontSize', 15);
zlabel('Height (m)','FontSize', 15);
axis([0 120 -30 30 0 30])
grid on


 
% Sub function to calculate x,y,vx,vz,mass of the system, the air volume, mass of air
function A = calcAccel(t,y,theta,l,v_wind)
c = getConst(); % Call the getConst function to get all constants


d = [y(1);y(2);y(3)]; % positions x and z
v = [y(4);y(5);y(6)]; % velocities x and z
v = v - v_wind;
m = y(7); % mass of the system 
V_air = y(8); % volume 
m_air = y(9); % mass of air
 
M_v = sqrt(v(1)^2 + v(2)^2 + v(3)^2); % magnitude of velocity 
h = v./M_v; % the heading vector (unit vector for x and z)
 
% Set initial conditions for the heading vector. 
if sqrt(d(1)^2 + d(2)^2 + d(3)^2) <= l
    h = [cosd(theta); 0; sind(theta)];
end
 
At = pi*c.D_thrust^(2)/4; % Cross sectional area of the throat
Ab = pi*c.D_bottle^(2)/4; % Cross sectional area of the bottle 
 
P_end = c.P0*(c.V_air0/c.Vbottle)^c.gamma; % Absolute pressure at the time all the water is expelled
P = P_end*(m_air/c.m_air0)^c.gamma; % Absolute pressure at any time in Phase2
    
% Phase 1: Water expulsion
if  V_air < c.Vbottle % The mass of the system should be greater than m_phase1
    Pa = c.P0*(c.V_air0/V_air)^c.gamma; % absolute pressure at any time t
    delta_P = Pa - c.P_atm; 
    V_e = sqrt(2*delta_P/c.rho_water); % pressure at the exit 
    
    mass = m; % the mass of the system
    fDrag = -0.5*c.rho_air*M_v^2*c.CD*Ab.*h; % drag forces (x and z) in Phase 1
    fThrust = 2*c.Cd*At*delta_P.*h; % thrust forces (x ans z)in Phase 1
    fGrav = [0; 0; -mass*c.g]; % gravitational force in Phase 1
    v = [y(4);y(5);y(6)]; % store the values of velocity at times in Phase 1
    m = -c.Cd*At*sqrt(2*c.rho_water*delta_P); % the rate of change in mass of the system with time
    V_air = c.Cd*At*V_e; % the rate of chang in voume of air with time
    m_air = 0; % no change in mass of air in this phase
 
% Phase 2: Water Exhasted 
elseif P > c.P_atm
    P_cri = P*(2/(c.gamma+1))^(c.gamma/(c.gamma-1)); % Critical pressure
    rho_air_inside = m_air/c.Vbottle; % since air is compressible the air density varies at any time in Phase2 
    T = P/(rho_air_inside*c.R); % Accordingly, the temperature also varies at any time
    
    if P_cri > c.P_atm % if critical pressure is greater than the atm pressure
        % Assume the Mach number at the exit is 1 (M_e = 1: sonic flow)
        T_e = T*(2/(c.gamma+1)); % the exit temp 
        P_e = P_cri; % the exit pressure should be the same as the critical pressure
        rho_e = P_e/(c.R*T_e); % the exit air density
        V_e = sqrt(c.gamma*c.R*T_e); % the exit velocity 
    else  
        M_e = ((2/(c.gamma-1))*((P/c.P_atm)^((c.gamma-1)/c.gamma) - 1))^(0.5); % the Mach number at the exit
        T_e = T/(1+((c.gamma-1)*M_e^(2)/2)); % the exit temp
        rho_e = c.P_atm/(c.R*T_e); % the exit density
        P_e = c.P_atm; % the exit pressure should be equal to the atm pressure
        V_e = M_e*sqrt(c.gamma*c.R*T_e); % the exit velocity
    end
    
    mass = m; 
    fGrav = [0;0; -mass*c.g]; % gravitational force
    fDrag = -0.5*c.rho_air*M_v^2*c.CD*Ab.*h; % drage foces in x and z
    fThrust = (c.Cd*rho_e*At*V_e^2 + (c.P_atm - P_e)*At).*h; % thrust foces in x and z
    v = [y(4);y(5);y(6)]; % current velocities in x and z
    m = -c.Cd*rho_e*At*V_e; % the rate of change in mass of the system
    V_air = 0; % No change in volume of air inside the bottle
    m_air = -c.Cd*At*rho_e*V_e;
    
% Phase 3: Ballistic
else
    mass = m; % the mass of the system should be equal to mass of the bottle
    fGrav = [0;0;-mass*c.g]; % Gravitational force
    fDrag = -0.5*c.rho_air*M_v^2*c.CD*Ab.*h; % drag forces in x and z
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
        m = 0.131;
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
c.CD = 0.1; % drag coefficient
c.rho_air = 0.961; % ambient air density (kg/m^3)
c.Vbottle = 0.002; % volume of empty bottle (m^3)
c.P_atm = 12.1*6894.76; % atmospheric pressure (Pa) (1 psi = 6894.76Pa)
c.gamma = 1.4; % ratio of the specific heat for air
c.rho_water = 997;  % density of wter (kg/m^3)
c.D_thrust = 0.021; % diameter of throat (m)
c.D_bottle = 0.105; % diameter of bottle (m)
c.R = 287; % air gas constant (J/kgK)
c.m_bottle = 0.129; % mass of empty 2-litter bottle (kg)
c.P_gage = 40*6894.76; % initial gage pressure (Pa) (1 psi = 6894.76Pa) 
c.m_water0 = 0.6; % initial mass of water (kg)
c.Vwater = c.m_water0/c.rho_water; % initial volume of water (kg)
c.T_air = 2 + 273.15; %K initial temperature (K)
c.P0 = c.P_atm + c.P_gage; % initial absolute pressuare (Pa)
c.V_air0 = c.Vbottle - c.Vwater; % initial volume of air inside the bottle (m^3)
c.m_air0 = c.P0*c.V_air0/(c.R*c.T_air); % initial mass of air (kg)
c.m0 = c.m_bottle + c.m_air0 + c.m_water0; % initial mass of the system (kg)

end
 

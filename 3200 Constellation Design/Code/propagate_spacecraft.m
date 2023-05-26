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
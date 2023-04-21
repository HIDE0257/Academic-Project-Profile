function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
a0_t = a0_t*180/pi;         % [/rad] Lift slope at the tips (NACA0012)
a0_r = a0_r*180/pi;         % [/rad] Lift slope at the root (NACA2412)
aero_t = aero_t*pi/180;     % [rad] Zero-lift angle of attack at the tips (NACA0012)
aero_r = aero_r*pi/180;     % [rad] Zero-lift angle of attack at the root (NACA2412)
geo_t = geo_t*pi/180;       % [rad] Geometric angle of attack at the tips (NACA0012)
geo_r = geo_r*pi/180;       % [rad] Geometric angle of attack at the root (NACA2412)
theta = @(x) x*pi/(2*N);    % [rad] Theta
alpha_geo = @(x) geo_r + (geo_t - geo_r)*cos(x);            % [rad] Varied Geometric anlge of attack
alpha_zerolift = @(x) aero_r + (aero_t - aero_r)*cos(x);    % [rad] Varied zero-lift angle of attack
c = @(x) c_r + (c_t - c_r)*cos(x);      % [ft] Varied chord line of wings
a0 = @(x) a0_r + (a0_t - a0_r)*cos(x);  % [/rad] Varied lift slope
s = b*(c_t + c_r)/2;        % [ft^2] Surface area of wings
AR = b^2/s;                 % Aspect ratio
 
for i = 1:N
    % F 1xN vector
    F_i(i,1) = alpha_geo(theta(i)) - alpha_zerolift(theta(i));
    for j = 1:N
        odd = 2*j-1;    % Odd terms
        % K NxN matrix 
        K_ij(i,j) = 4*b*sin(odd*theta(i))/(a0(theta(i))*c(theta(i))) + (odd*sin(odd*theta(i))/sin(theta(i)));
    end
end

% Coefficient 1xN vector A
A_j = K_ij\F_i;

% Odd 1xN-1 vector that contains only odd terms starting 3
for k = 2:N
    odd(k-1) = 2*k-1;
end
e = (1 + sum(odd(:).*(A_j(2:N)./A_j(1)).^2))^(-1);  % Efficiency factor 
c_L = A_j(1)*pi*AR;         % Lift coefficient 
c_Di = c_L^2/(2*pi*e*AR);   % Induced drag coefficient 
end
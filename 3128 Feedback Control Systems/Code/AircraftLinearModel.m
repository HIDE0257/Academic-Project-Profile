function [Alon, Blon, Alat, Blat] = AircraftLinearModel(trim_definition, trim_variables, aircraft_parameters)
%
% state = [x^E, y^E, Z^E, roll, pitch, yaw, u^E, v^E, w^E, p, q, r]';
%
%


u0 = trim_definition(1);
h0 = trim_definition(2);

alpha0 = trim_variables(1);
de0 = trim_variables(2);
dt0 = trim_variables(3);

theta0 = alpha0;

ap = aircraft_parameters;
rho = stdatmo(h0);



%%%%%%%%%%%%%%%%%%%%%%
%%% Longitudinal
%%%%%%%%%%%%%%%%%%%%%%
%%%% Trim values
CW0 = ap.W/((1/2)*rho*u0^2*ap.S);
CL0 = CW0*cos(theta0);
CD0 = ap.CDmin + ap.K*(CL0-ap.CLmin)^2;
CT0 = CD0 + CW0*sin(theta0);
%%%% Nondimensional stability derivatives in body coordinates
dTdu = dt0*ap.Cprop*ap.Sprop*(ap.kmotor-2*u0+dt0*(-2*ap.kmotor+2*u0));
CXu = dTdu/(.5*rho*u0*ap.S)-2*CT0;
CDu = 0;
CLu = 0;
Cmu = 0;
CZu = -CLu;
CZalpha = -CD0 - ap.CLalpha;
CXalpha = CL0*(1-2*ap.K*ap.CLalpha);
CZalphadot = -ap.CLalphadot;
CZq = -ap.CLq;
% Longitudinal dimensional stability derivatives (from Etkin and Reid)
Xu = rho*u0*ap.S*(CW0*sin(theta0) + .5*CXu);
Zu = rho*u0*ap.S*(-CW0*cos(theta0) + .5*CZu);
Mu = (1/2)*rho*u0*ap.c*ap.S*Cmu;
Xw = (1/2)*rho*u0*ap.S*CXalpha;
Zw = (1/2)*rho*u0*ap.S*CZalpha;
Mw = (1/2)*rho*u0*ap.c*ap.S*ap.Cmalpha;
Xq = 0;
Zq = (1/4)*rho*u0*ap.c*ap.S*CZq;
Mq = (1/4)*rho*u0*ap.c^2*ap.S*ap.Cmq;
Xwdot = 0;
Zwdot = (1/4)*rho*ap.c*ap.S*CZalphadot;
Mwdot = (1/4)*rho*ap.c^2*ap.S*ap.Cmalphadot;
CDde = 2*ap.K*(CL0-ap.CLmin)*ap.CLde;
Cxde = -CDde*cos(alpha0)+ap.CLde*sin(alpha0);
Czde = -CDde*sin(alpha0)-ap.CLde*cos(alpha0);
% Elevator
Xde = 0.5*rho*(u0^2)*ap.S*Cxde;
Zde = 0.5*rho*(u0^2)*ap.S*Czde;
Mde = 0.5*rho*(u0^2)*ap.S*ap.c*ap.Cmde;
% Tail
CXdt = 2*dt0*((ap.Sprop*ap.Cprop)/ap.S)*(ap.kmotor/u0)^2;
CZdt = 0; % Thrust Alligned
Cmdt = 0; % Thrust Alligned
Xdt = 0.5*rho*(u0^2)*ap.S*CXdt;
Mdt = 0.5*rho*(u0^2)*ap.S*ap.c*Cmdt;
Zdt = 0.5*rho*(u0^2)*ap.S*CZdt;
% Matrices
% x_lat = [du, dw, dq, dtheta, dx, dz]
% u_lat = [delta_a, delta_r]
Alon1 = [  Xu/ap.m,           Xw/ap.m,               0,                          -ap.g*cos(theta0);
       Zu/(ap.m-Zwdot),   Zw/(ap.m-Zwdot),       (Zq + ap.m*u0)/(ap.m-Zwdot),      -(ap.m*ap.g*sin(theta0))/(ap.m-Zwdot);
       (Mu + (Mwdot*Zu)/(ap.m-Zwdot))/ap.Iy,         (Mw + (Mwdot*Zw)/(ap.m-Zwdot))/ap.Iy, ...
               (Mq + (Mwdot*(Zq+ap.m*u0))/(ap.m-Zwdot))/ap.Iy,      ((Mwdot*(-ap.m*ap.g*sin(theta0)))/(ap.m-Zwdot))/ap.Iy;
       0,              0,                  1,                             0];
Blon = [Xde/ap.m,    Xdt/ap.m;
       Zde/(ap.m -Zwdot), Zdt/(ap.m - Zwdot);
       (Mde/ap.Iy) + (Mwdot/ap.Iy)*(Zde/(ap.m-Zwdot)) , (Mdt/ap.Iy) + (Mwdot/ap.Iy)*(Zdt/(ap.m-Zwdot));
       0,  0;
       0,  0;
       0,  0];
Alon = [Alon1 zeros(4,2);...
       1 0 0 0 0 0;...
       0 1 0 -u0 0 0];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lateral
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lateral-directional dimensional stability derivatives
Yv = (1/2)*rho*u0*ap.S*ap.CYbeta;
Yp = (1/4)*rho*u0*ap.b*ap.S*ap.CYp;
Yr = (1/4)*rho*u0*ap.b*ap.S*ap.CYr;
Ydr = (1/2)*rho*u0^2*ap.S*ap.CYdr;
Yda = (1/2)*rho*u0^2*ap.S*ap.CYda;

Lv = (1/2)*rho*u0*ap.b*ap.S*ap.Clbeta;
Lp = (1/4)*rho*u0*ap.b^2*ap.S*ap.Clp;
Lr = (1/4)*rho*u0*ap.b^2*ap.S*ap.Clr;
Lda = (1/2)*rho*u0^2*ap.b*ap.S*ap.Clda;
Ldr = (1/2)*rho*u0^2*ap.b*ap.S*ap.Cldr;

Nv = (1/2)*rho*u0*ap.b*ap.S*ap.Cnbeta;
Np = (1/4)*rho*u0*ap.b^2*ap.S*ap.Cnp;
Nr = (1/4)*rho*u0*ap.b^2*ap.S*ap.Cnr;
Nda = (1/2)*rho*u0^2*ap.b*ap.S*ap.Cnda;
Ndr = (1/2)*rho*u0^2*ap.b*ap.S*ap.Cndr;

G = ap.Ix*ap.Iz-ap.Ixz^2;

G3=ap.Iz/G;
G4=ap.Ixz/G;
G8=ap.Ix/G;


% x_lat = [dv, dp, dr, dphi, dpsi, dy]
% u_lat = [delta_a, delta_r]
Alat1 = [Yv/ap.m, Yp/ap.m, (Yr/ap.m - u0), ap.g*cos(theta0);
    G3*Lv + G4*Nv, G3*Lp + G4*Np, G3*Lr + G4*Nr, 0;
    G4*Lv + G8*Nv, G4*Lp + G8*Np, G4*Lr + G8*Nr, 0;
    0, 1, tan(theta0), 0];

Blat = [Yda/ap.m Ydr/ap.m; G3*Lda+G4*Nda G3*Ldr+G4*Ndr; G4*Lda+G8*Nda G4*Ldr+G8*Ndr; 0 0; 0 0; 0 0];


Alat = [Alat1 zeros(4,2);...
        0 0 sec(theta0) 0 0 0;...
        1 0 0 0 u0*cos(theta0) 0];


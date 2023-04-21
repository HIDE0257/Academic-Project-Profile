function [Cl,Cd_w] = LinearizedTheory(c,alpha,M,e1,e2)
% c = chord length (m)
% alpha[vector] = angle of attack (deg)
% M = Mach number
% e1 = geometric angle of diamond airfoil in front
% e2 = geometric angle of diamond airfoil in back

    alpha = deg2rad(alpha);     % [rad] Angle of attack
    e1 = deg2rad(e1);           % [rad] epsilon1
    e2 = deg2rad(e2);           % [rad] epsilon2
    c1 = tan(e2)*c/(tan(e1) + tan(e2)); % [m] Chord length relative to epsilon1
    c2 = tan(e1)*c/(tan(e1) + tan(e2)); % [m] Chord length relative to epsilon2
    f1_u = @(x) tan(e1)^2*x;    % Integrated (dy_u/dx)^2 relative to epsilon1
    f2_u = @(x) tan(e2)^2*x;    % Integrated (dy_u/dx)^2 relative to epsilon2
    f1_l = @(x) tan(e1)^2*x;    % Integrated (dy_l/dx)^2 relative to epsilon1
    f2_l = @(x) tan(e2)^2*x;    % Integrated (dy_l/dx)^2 relative to epsilon2
    
    g_u = (1/c)*(f1_u(c1) + f2_u(c) - f2_u(c1));   
    g_l = (1/c)*(f1_l(c1) + f2_l(c) - f2_l(c1));
    Cl = 4.*alpha(:)./(sqrt(M.^2 - 1));     % Linearized lift coefficient 
    Cd_w = (2./sqrt(M.^2 - 1)).*(2.*alpha(:).^2 + g_u + g_l);   % Linearized wave drag coeff
end
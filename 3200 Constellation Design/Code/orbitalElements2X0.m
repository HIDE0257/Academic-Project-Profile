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
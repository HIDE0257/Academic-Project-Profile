% Function with constants for the two-element model 
function Constants2 = getConstants2()
Constants2.E = 10175000 ; % Young's Modulus [lb/in^2]
Constants2.w = 1; % Width of all members [in] 
Constants2.h = 1/8; % Thickness of fuselage member [in]
Constants2.A = Constants2.w*Constants2.h; % Area [in^2]
Constants2.L = 12; % End-shaker-to-start-tail span [in]
Constants2.rho = 0.0002505; % Density of bar materials [lb-sec^2/in^4]
Constants2.Mt = 1.131*Constants2.rho; % Mass of tail assembly [slugs]
Constants2.St = 0.5655*Constants2.rho; % First mass-moment of tail assembly wrt to B' 
Constants2.It = 23.124*Constants2.rho; % Second mass-moment of tail assembly wrt B'
Constants2.cM2 = (Constants2.rho*Constants2.A*Constants2.L)/(100800); % Coefficient for mass matrix 
Constants2.Izz = (Constants2.w*Constants2.h^3)/12; % Inertia about z-axis
Constants2.cK2 = (4*Constants2.E*Constants2.Izz)/(Constants2.L^3); % Coefficient for stiffness matrix
Constants2.cM4 = (Constants2.rho*Constants2.A*Constants2.L)/806400; % Coefficient for four element master mass equation 
Constants2.cK4 = (8*Constants2.E*Constants2.Izz)/(Constants2.L^3); % Coefficient for four element master stiffness equation
end
function zerolift_TAT = Thin_Airfoil_Theory(m,p)
% NACA(mptt)
m = m./100;     % Maximum cambered line in 100th 
p = p./10;      % Location of maximum cambered line in 10th
limits2412 = [0, acos(1/5); acos(1/5), pi]; % Boundaries for integration 
limits4412 = [0, acos(1/5); acos(1/5), pi]; % Boundaries for integration

% Functions of inside the integration for T.A.T zerolift angle of attack 
F1_2412 = @(x) (2.*m(2)./p(2)).*(1 - (1 - cos(x))./(2.*p(2))).*(cos(x) - 1);
F2_2412 = @(x) (2.*m(2)./(1 - p(2)).^2).*(p(2) - (1 - cos(x))./2).*(cos(x) - 1);
F1_4412 = @(x) (2.*m(3)./p(3)).*(1 - (1 - cos(x))./(2.*p(3))).*(cos(x) - 1);
F2_4412 = @(x) (2.*m(3)./(1 - p(3)).^2).*(p(3) - (1 - cos(x))./2).*(cos(x) - 1);

% Zerolift angle of attack of NACA0012 (Symmetric)
zerolift_TAT1 = 0; 
% Zerolift angle of attack of NACA2412
zerolift_TAT2 = -(1/pi).*(integral(F1_2412,limits2412(1,1),limits2412(1,2)) + integral(F2_2412,limits2412(2,1),limits2412(2,2)));
% Zerolift angle of attack of NACA4412
zerolift_TAT3 = -(1/pi).*(integral(F1_4412,limits4412(1,1),limits4412(1,2)) + integral(F2_4412,limits4412(2,1),limits4412(2,2)));
% Combined them above into a vector form 
zerolift_TAT = [zerolift_TAT1, zerolift_TAT2, zerolift_TAT3].*180/pi;
end
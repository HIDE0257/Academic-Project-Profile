%% DCM_EtoB
function DCM_EtoB = T_EtoB(i,Omega,w)
% i = [rad] inclination
% Omega = [rad] right ascension
% w = [rad] argument of periapsis 

R1 = [cos(Omega) sin(Omega) 0; -sin(Omega) cos(Omega) 0; 0 0 1]; 
R2 = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)]; 
R3 = [cos(w) sin(w) 0; -sin(w) cos(w) 0; 0 0 1]; 

DCM_EtoB = [R3*R2*R1];

end
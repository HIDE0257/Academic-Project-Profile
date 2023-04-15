% Assumptions *********************************
% - SLUF ( L = W )
% - STD ATM @ 1500m 
% 

% Requirements ********************************
R_req = [125,175];      % (m) Min - Max of glider range
V_Rreq = [7, 12];       % (m/s) minimum velocity for range
E_req = [13, 20];       % (s) maximum time of flight
theta_req = [10, 20];   % (deg) Elevetor pitch control
vol_HT = [0.3,0.6];     % Volume coefficient of horizontal tail
vol_VT = [0.02,0.05];   % Volume coefficient of vertical tail
SM_req = [0.15,0.25];       % Longitudinal Stability

% Constants ***********************************
chord = 0.18;       % (m) Depth of the wing
%S_wet = 0.49;       %((m^2) (2*S_w) + (2*S_ht) + (2*S_vt) + S_f;
b = 1;            % (m) Length of wing span
S_ht = 0.036;        % (m^2) Surface area of horizontal tail
S_vt = 2*0.007;        % (m^2) Surface area of vertical tail
S_ref = b*chord;    % (m^2) Reference area
S_f = 0.033;        % (m^2) Fuselage area 
S_wet = 2*S_ref + 2*S_ht + 2*S_vt + S_f;
AR = (b^2)/S_ref;

% e0 = 1.78*(1 - (0.045*(AR^0.68))) - 0.64;           % Efficiency factor
e0 = 0.7;
h = 17.5;           % (m) Height where a glider is thrown at 
rho = 1.0581;       % (kg/m^3) At an altitude of 1500m
C_fe = 0.003;       % Equivalent coefficient of skin friction
ws = 0.295*9.81;    % (N/m^2) Dollar Tree Foam's approximate specific weight
W_cam = 0.16*9.81;  % (N) Weight of camera
CL_max = 0.8;       % Maximum lift coefficient from the data file

W_TO = ws*S_wet + W_cam; %(S_f+S_w+S_ht+S_vt)
k = 1/(AR*pi*e0);
wing_load = W_TO/S_ref;

%**********************************************
% CD & CL
CD0 = C_fe*(S_wet/S_ref);       % Paracite drag
CD_R = 2*CD0;                 % Minimum drag polar
CL_R = sqrt(CD0/k);           % Maximum CL : CD0 = kCL^2
ratio_LD_R = CL_R/CD_R;       % Ratio of lift to drag 
V_stall = sqrt(2*wing_load/(CL_max*rho)); % Stall velocity

% Range & Velocity at the R_max
R = h*ratio_LD_R;
V_R = sqrt(2*wing_load/(rho*CL_R));   % Velocity for the maximum range

% Endurance & Velocity at the E_min
CL_E = sqrt(3*CD0/k);         % Minimim CL to get minimum sink velocity: 3CD0 = kCL^2
CD_E = 4*CD0;
ratio_LD_E = CL_R/CD_E;
theta_E = (1/ratio_LD_R);
V_sink = V_R*sin(theta_E);
V_E = sqrt(2*wing_load/(rho*CL_E));
E = h/V_sink;

theta_deg = theta_E*180/pi;
% 
% for t = 0:E
%     x = V_R*t;
%     y = -V_sink*t + h;
%     if y >= 0
%         figure(5)
%         hold on
%         plot(x,y,'ro');
%         title('Flight Simulation Beased on Our Actual Design');
%         xlabel('Displacement (m)');
%         ylabel('Height (m)');
%         axis([0 200 0 20]);
%         drawnow
%     else 
%         break;
%     end
%     
% end

% Display key values 
fprintf('                Range:   %.1f m\n', R);
fprintf('             Ref Area:   %.2f m^2\n', S_ref);
fprintf('          Wetted Area:   %.2f m^2\n', S_wet);
fprintf('         Wing Loading:   %.1f\n', wing_load);
fprintf('         Aspect Ratio:   %.1f\n', AR);
fprintf('    L/D max for Range:   %.1f\n', ratio_LD_R);
fprintf('L/D max for Endurance:   %.1f\n', ratio_LD_E);
fprintf('   Velocity for Range:   %.1f m/s\n', V_R);
fprintf('            Endurance:   %.1f s\n', E);
fprintf('       Elevator Pitch:   %.1f deg\n', theta_deg);

clear all;
clc
% Assumptions *********************************
% - SLUF ( L = W )
% - STD ATM @ 1500m

% Requirements ********************************
R_req = [125,175];      % (m) Min - Max of glider range
V_Rreq = [7, 12];       % (m/s) minimum velocity for range
E_req = [13, 25];       % (s) maximum time of flight
theta_req = [10, 20];   % (deg) Elevetor pitch control
vol_HT_req = [0.3,0.6];     % Volume coefficient of horizontal tail
vol_VT_req = [0.02,0.05];   % Volume coefficient of vertical tail

% Constants ***********************************
e0 = 0.7;           % Efficiency factor
h = 17.5;           % (m) Height where a glider is thrown at
rho = 1.0581;       % (kg/m^3) At an altitude of 1500m
b = 1;              % (m) Length of wing span
C_fe = 0.003;       % Equivalent coefficient of skin friction
ws = 0.295*9.81;    % (N/m^2) Dollar Tree Foam's approximate specific weight
W_pay = 0.16*9.81;  % (N) Weigth of camera
CL_max = 0.8;       % Maximum lift coefficient from the data file
S_ht = 0.036;        % (m^2) Surface area of horizontal tail
S_vt = 0.014;        % (m^2) Surface area of vertical tail
S_f = 0.033;        % (m^2) Fuselage area 

chord = linspace(0.1,1,100);   % (m) Chord length

for i = 1:length(chord)
        S_ref(i) = b*chord(i);     % (m^2) Refference area (Wing body area)
        S_wet(i) = 2*S_ref(i) + 2*S_ht + 2*S_vt + S_f;
        W_TO(i) = ws*S_wet(i) + W_pay; % (N) Estimated weigth of the aircraft
        AR(i) = (b^2)/S_ref(i);     % Aspect ratio
        k(i) = 1/(AR(i)*pi*e0);     % Coefficient for L and D
        wing_load(i) = W_TO(i)/S_ref(i);        % (N/m^2) Wing loading
        
        %**********************************************
        % CD & CL
        CD0(i) = C_fe*(S_wet(i)/S_ref(i));        % Paracite drag
        CD_R(i) = 2*CD0(i);                     % Minimum drag polar for range
        CL_R(i) = sqrt(CD0(i)/k(i));            % CL for range : CD0 = kCL^2
        ratio_LD_R(i) = CL_R(i)/CD_R(i);     % Ratio of lift to drag for range (L/D_max)
        V_stall(i) = sqrt(2*wing_load(i)/(CL_max*rho)); % Limit velocity
        
        
        % Range & Velocity at the R_max
        R(i) = h*ratio_LD_R(i);     % (m) Range
        V_R(i) = sqrt(2*wing_load(i)/(rho*CL_R(i)));   % (m/s) Velocity for range
        % It's assumed that V_R is only in x direction since the SLUP assumption
        
        % Endurance & Velocity at the E_min
        CD_E(i) = 4*CD0(i);                     % Minimum drag polar for endurance
        CL_E(i) = sqrt(3*CD0(i)/k(i));          % CL for endurance: 3CD0 = kCL^2
        ratio_LD_E(i) = CL_R(i)/CD_E(i);      % Ratio of lift to drag for endurance (L/D_max)
        theta_R(i) = 1/ratio_LD_R(i);         % (rad) Sink Angle
        theta_E(i) = 1/ratio_LD_E(i);
        V_sink(i) = V_R(i)*sin(theta_R(i));   % (m/s) Sink velocity
        E(i) = h/V_sink(i);                     % (s) Endurance
        V_E(i) = sqrt(2*wing_load(i)/(rho*CL_E(i)));  % (m/s^2) Velocity for endurance
        
        theta_deg(i) = theta_E(i)*180/pi;         % (deg) Sink Angle
        
        if R_req(1) <= R(i) && R(i) <= R_req(2)
%             if V_Rreq(1) <= V_R(i) && V_R(i) <= V_Rreq(2)
                if E_req(1) <= E(i) && E(i) <= E_req(2)
%                      if 7 <= theta_deg(i) && theta_deg(i) <= theta_req(2)
                        if 7.14 <= ratio_LD_R(i) && ratio_LD_R(i) <= 10
                            best_R(i) = R(i);
                            best_wl(i) = wing_load(i);
                        end
%                      end
                end
%             end
        end
    
end

c = 10;

best_chord = chord(c);
best_Sref = S_ref(c);
best_Swet = S_wet(c);

fprintf('Wing Reference Area: %.2f m^2\n', best_Sref);
fprintf('        Wetted Area: %.2f m^2\n', best_Swet); 


figure(1)       % Plot of range & required CL vs wing loading
hold on
subplot(2,1,1);
yyaxis left
plot(wing_load(:),R(:));
yline([R_req(1),R_req(2)],'--',{'Min','Max'});
xline(wing_load(c),'r', {'Desired Wing Loading'});
title('Actual Design: Wing Loading and CL Required vs Range');
xlabel('Wing Loading W/S (N/m^2)');
ylabel('Range (m)');
axis([0 45 160 450]);

yyaxis right
plot(wing_load(:),CL_R(:));
ylabel('Lift Coefficient Required');
yline(CL_max,'--g',{'Max CL'});
axis([0 45 0.1 0.8]);
grid on

subplot(2,1,2);
yyaxis left
plot(wing_load(:),E(:));
yline([E_req(1),E_req(2)],'--',{'Min','Max'});
xline(wing_load(c),'r', {'Desired Wing Loading'});
title('Actual Design: Wing Loading and CL Required vs Endurance');
xlabel('Wing Loading W/S (N/m^2)');
ylabel('Endurance (s)');
axis([0 45 12 45]);

yyaxis right
plot(wing_load(:),CL_E(:));
ylabel('Lift Coefficient Required');
yline(CL_max,'--g',{'Max CL'});
axis([0 45 0.2 0.9]);
grid on
hold off

figure(2)       % Plot of endurance & required CL vs wing loading
hold on
plot(wing_load(:),V_R(:));
plot(wing_load(:),V_E(:));
plot(wing_load(:),V_stall(:),'--');
xline(wing_load(c),'r', {'Desired Wing Loading'});
yline([V_Rreq(1),V_Rreq(2)],'--',{'Min','Max'});
title('Actual Design: Velocity Required for Max R & E');
xlabel('Wing Loading W/S (N/m^2)');
ylabel('Velocity Required (m/s)');
legend('Velocity to achieve Rmax','Velocity to achieve Emax','Limit Velocity for Maximum CL');
axis([0 40 5 15]);
grid on
hold off


figure(3)
hold on
xline(wing_load(c),'r', {'Desired Wing Loading'});
plot(wing_load(:),W_TO(:));
title('Actual Design: Variation of Aircraft Weight with Wing Loading');
ylabel('Aircraft Weight (N)');
xlabel('Wing Loading (N/m^2)');
axis([0 20 2.5 5]);
hold off


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
C_fe = 0.003;      % Equivalent coefficient of skin friction
ws = 0.295*9.81;    % (N/m^2) Dollar Tree Foam's approximate specific weight
W_pay = 0.16*9.81;  % (N) Weigth of camera
CL_max = 0.8;       % Maximum lift coefficient from the data file
S_ht = 0.04;        % (m^2) Surface area of horizontal tail
S_vt = 0.011;       % (m^2) Surface area of vertical tail


chord = linspace(0.01,1,100);   % (m) Chord length
S_wet = linspace(0.07,10,100);  % (m^2) Wetted area

for i = 1:length(chord)
    for j = 1:length(S_wet)
        S_ref(i) = b*chord(i);     % (m^2) Refference area (Wing body area)
        W_TO(j) = ws*S_wet(j) + W_pay; % (N) Estimated weigth of the aircraft
        AR(i) = (b^2)/S_ref(i);     % Aspect ratio
        k(i) = 1/(AR(i)*pi*e0);     % Coefficient for L and D
        wing_load(i,j) = W_TO(j)/S_ref(i);        % (N/m^2) Wing loading
        
        %**********************************************
        % CD & CL
        CD0(i,j) = C_fe*(S_wet(j)/S_ref(i));        % Paracite drag
        CD_R(i,j) = 2*CD0(i,j);                     % Minimum drag polar for range
        CL_R(i,j) = sqrt(CD0(i,j)/k(i));            % CL for range : CD0 = kCL^2
        ratio_LD_R(i,j) = CL_R(i,j)/CD_R(i,j);     % Ratio of lift to drag for range (L/D_max)
        V_stall(i,j) = sqrt(2*wing_load(i,j)/(CL_max*rho)); % Limit velocity
        
        
        % Range & Velocity at the R_max
        R(i,j) = h*ratio_LD_R(i,j);     % (m) Range
        V_R(i,j) = sqrt(2*wing_load(i,j)/(rho*CL_R(i,j)));   % (m/s) Velocity for range
        % It's assumed that V_R is only in x direction since the SLUP assumption
        
        % Endurance & Velocity at the E_min
        CD_E(i,j) = 4*CD0(i,j);                     % Minimum drag polar for endurance
        CL_E(i,j) = sqrt(3*CD0(i,j)/k(i));          % CL for endurance: 3CD0 = kCL^2
        ratio_LD_E(i,j) = CL_R(i,j)/CD_E(i,j);      % Ratio of lift to drag for endurance (L/D_max)
        theta_R(i,j) = 1/ratio_LD_R(i,j);         % (rad) Sink Angle
        theta_E(i,j) = 1/ratio_LD_E(i,j);
        V_sink(i,j) = V_R(i,j)*sin(theta_R(i,j));   % (m/s) Sink velocity
        E(i,j) = h/V_sink(i,j);                     % (s) Endurance
        V_E(i,j) = sqrt(2*wing_load(i,j)/(rho*CL_E(i,j)));  % (m/s^2) Velocity for endurance
        
        theta_deg(i,j) = theta_R(i,j)*180/pi;         % (deg) Sink Angle
        
        if R_req(1) <= R(i,j) && R(i,j) <= R_req(2)
            if V_Rreq(1) <= V_R(i,j) && V_R(i,j) <= V_Rreq(2)
                if E_req(1) <= E(i,j) && E(i,j) <= E_req(2)
%                     if 7 <= theta_deg(i,j) && theta_deg(i,j) <= theta_req(2)
                        if 7.14 <= ratio_LD_R(i,j) && ratio_LD_R(i,j) <= 10
                            best_R(i,j) = R(i,j);
                            best_wl(i,j) = wing_load(i,j);
                        end
%                     end
                end
            end
        end
    end
end

[m,n] = find(best_wl(:,:));

r = 18;
c = 5;

best_Sref = S_ref(r);
best_Swet = S_wet(c);

fprintf('Wing Reference Area: %.2f m^2\n', best_Sref);
fprintf('        Wetted Area: %.2f m^2\n', best_Swet); 
row = r;
col = c;

figure(1)       % Plot of range & required CL vs wing loading
hold on
subplot(2,1,1);
yyaxis left
plot(wing_load(row,:),R(row,:));
yline([R_req(1),R_req(2)],'--',{'Min','Max'});
xline(wing_load(row,col),'r', {'Desired Wing Loading'});
title('Actual Design: Wing Loading and CL Required vs Range');
xlabel('Wing Loading W/S (N/m^2)');
ylabel('Range (m)');
axis([0 45 120 450]);

yyaxis right
plot(wing_load(row,:),CL_R(row,:));
ylabel('Lift Coefficient Required');
yline(CL_max,'--g',{'Max CL'});
axis([0 45 0.1 0.8]);
grid on

subplot(2,1,2);
yyaxis left
plot(wing_load(row,:),E(row,:));
yline([E_req(1),E_req(2)],'--',{'Min','Max'});
xline(wing_load(row,col),'r', {'Desired Wing Loading'});
title('Actual Design: Wing Loading and CL Required vs Endurance');
xlabel('Wing Loading W/S (N/m^2)');
ylabel('Endurance (s)');
axis([0 45 10 45]);

yyaxis right
plot(wing_load(row,:),CL_E(row,:));
ylabel('Lift Coefficient Required');
yline(CL_max,'--g',{'Max CL'});
axis([0 45 0.2 0.9]);
grid on
hold off

figure(2)       % Plot of endurance & required CL vs wing loading
hold on
plot(wing_load(row,:),V_R(row,:));
plot(wing_load(row,:),V_E(row,:));
plot(wing_load(row,:),V_stall(row,:),'--');
xline(wing_load(row,col),'r', {'Desired Wing Loading'});
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
xline(wing_load(row,col),'r', {'Desired Wing Loading'});
plot(wing_load(:,col),W_TO(:));
title('Actual Design: Variation of Aircraft Weight with Wing Loading');
ylabel('Aircraft Weight (N)');
xlabel('Wing Loading (N/m^2)');
axis([0 45 2 30]);
hold off

figure(4)
hold on
xline(wing_load(row,col),'r', {'Desired Wing Loading'});
plot(wing_load(:,col),W_TO(:).*100);
title('Actual Design: Variation of Cost with Wing Loading ($100 /N)');
ylabel('Aircraft Weight (N)');
xlabel('Aircraft Cost ($)');
axis([0 45 200 3000]);
hold off

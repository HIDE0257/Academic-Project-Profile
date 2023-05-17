clear all;
close all;

recuv_tempest;

Va_trim = 21;
h_trim = 1800;

wind_inertial = [0;0;0];

trim_definition = [Va_trim; h_trim];
tfinal = 200;

%%% Use full minimization to determine trim
[trim_variables, fval] = CalculateTrimVariables(trim_definition, aircraft_parameters);
[trim_state, trim_input]= TrimStateAndInput(trim_variables, trim_definition);
[Alon, Blon, Alat, Blat] = AircraftLinearModel(trim_definition, trim_variables, aircraft_parameters);

[Evec_lat,Eval_lat] = eig(Alat);
Evec_roll = Evec_lat(:,3); 
Evec_roll = real(Evec_roll./Evec_roll(4))*deg2rad(5);
init_Roll12 = [0;Evec_roll(6);0;Evec_roll(4);0;Evec_roll(5);0;Evec_roll(1);0;Evec_roll(2);0;Evec_roll(3)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper code for Problem 3. Look inside the AircraftEOMControl function
% for hints on how to set up controllers for Problem 2 and Problem 4.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% K values 
Kr = -0.3; 
Ka = -5.5; % Between -4.5 and -5.5(better) is good
Kp = 3;
Kq = -0.15;
Kth = -10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Control Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi_c = deg2rad(5);
theta_c = deg2rad(5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yaw Damper
% x_lat = [v, p, r, phi, psi, y] 
[num_yaw, den_yaw] = ss2tf(Alat(1:4,1:4), Blat(1:4,2), [0 0 -1 0], 0); 

num_c_yaw = [0 0 Kr 0];
den_c_yaw = 1; 

yaw_cl = feedback(tf(conv(num_c_yaw, num_yaw), conv(den_c_yaw, den_yaw)),1);
[num_cl_yaw, den_cl_yaw] = tfdata(yaw_cl,'v');
roots(den_cl_yaw);
Nonlinear_yaw_ss = ss(Alat(1:4,1:4), Blat(1:4, 2), [0 0 1 0], 0);

figure(100)
subplot(1,2,1)
rlocus(num_yaw, den_yaw)
sgrid
title('Root Locus Yaw Damper')

subplot(1,2,2)
step(yaw_cl,50); hold on
step(Nonlinear_yaw_ss,50); hold on 
title('Roll Controller (Roll Rate to Aileon)')
legend(sprintf('Controlled Model (Kr=%.2f)', Kr), 'No-Control Model')

K_yaw = [0 0 Kr 0 0 0];
New_Alat_yaw = Alat - Blat(:,2)*K_yaw;
[New_Evec_yaw,New_Eval_yaw] = eig(New_Alat_yaw);

% Non-Control and Control Simulations
Evec_Yaw = New_Evec_yaw(:,4);
Evec_Yaw = real(Evec_Yaw/Evec_Yaw(2))*deg2rad(5); 
init_Yaw12_c = [0;Evec_Yaw(6);0;Evec_Yaw(4);0;Evec_Yaw(5);0;Evec_Yaw(1);0;Evec_Yaw(2);0;Evec_Yaw(3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Full sim in ode45
aircraft_state0 = trim_state + init_Yaw12_c;
aircraft_state0_y_c= trim_state + init_Yaw12_c; 
control_input0 = trim_input; % [de, da, dr, dt];
TSPAN = [0 200];

[T_noControl_y,Y_noControl_y] = ode45(@(t,y) AircraftEOM(t,y,control_input0,wind_inertial,aircraft_parameters),TSPAN,aircraft_state0);
[T_Control_y,Y_Control_y] = ode45(@(t,y) AircraftEOMControl(t,y,control_input0,wind_inertial,aircraft_parameters,0,0,0,0,Kr),TSPAN,aircraft_state0_y_c);
% RollControl(ka, pc, dp) = da    pc = roll rate command, dp = delta roll
% rate

for i = 1:length(T_noControl_y)
    Uout_noControl_y(i,:) = control_input0';
end
for i = 1:length(T_Control_y)
    dr_nonlinear_y(i,1) = Kr*Y_Control_y(i,12); % ka, kp,phi_c,phi,p)
    Uout_nonlinear_roll_y(i,:) = control_input0' + [0 0 dr_nonlinear_y(i) 0];
end

%% 2) Pitch Angle Control
%%% Transfer function from elevator to pitch angle
% x_lon = [u, w, q, theta, x, z]'   K_lon = [k_u, k_w, k_q, k_theta]'
% k_theta = proportional gain 
% k_q = derivative gain
[num_elev2pitch, den_elev2pitch] = ss2tf(Alon(1:4,1:4), Blon(1:4,1), [0 0 0 1],0);

%%% Controller
num_c = [0 0 Kq Kth];
den_c = 1;

K = [0 0 Kq Kth 0 0];

%%% Closed loop transfer function
pitch_cl = feedback(tf(conv(num_c, num_elev2pitch), conv(den_c, den_elev2pitch)), 1);
[num_cl, den_cl] = tfdata(pitch_cl,'v');
New = Alon - Blon(1:6, 1)*K;
[New_eigenVec, New_eigenVal] = eig(New);
[New_freq, New_damp] = damp(New);
[originalVec, originalVal] = eig(Alon);
[original_frq, original_damp] = damp(Alon);
openLoop = ss(Alon, Blon(:,1), [0 0 0 1 0 0], 0);

Evec_Pitch_c = New_eigenVec(:,5); 
Evec_Pitch_c = real(Evec_Pitch_c/Evec_Pitch_c(1))*deg2rad(5);
init_Pitch12_c = [Evec_Pitch_c(5); 0; Evec_Pitch_c(6); 0;Evec_Pitch_c(4);0;Evec_Pitch_c(1); 0;Evec_Pitch_c(2);0;Evec_Pitch_c(3);0];

Evec_Pitch = originalVec(:,5); 
Evec_Pitch = real(Evec_Pitch./Evec_Pitch(1))*deg2rad(5);
init_Pitch12 = [Evec_Pitch(5); 0; Evec_Pitch(6); 0;Evec_Pitch(4);0;Evec_Pitch(1); 0;Evec_Pitch(2);0;Evec_Pitch(3);0];

% Poles of the closed loop (linear system)
figure (101)
hold on;
    step(pitch_cl, 200);
    step(openLoop, 200);
    legend('Linear Model', 'Nonlinear Model')
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Full sim in ode45
aircraft_state0 = trim_state + init_Pitch12;
aircraft_state0_pitch_c= trim_state + init_Pitch12_c; 
control_input0 = trim_input; % [de, da, dr, dt];
TSPAN = [0 200];

[T_noControl_pitch,Y_noControl_pitch] = ode45(@(t,y) AircraftEOM(t,y,control_input0,wind_inertial,aircraft_parameters),TSPAN,aircraft_state0);
[T_Control_pitch,Y_Control_pitch] = ode45(@(t,y) AircraftEOMControl(t,y,control_input0,wind_inertial,aircraft_parameters,Kth,Kq,0,0,0),TSPAN,aircraft_state0_pitch_c);
% RollControl(ka, pc, dp) = da    pc = roll rate command, dp = delta roll
% rate

% control_input = [de, da, dr, dt];
for i = 1:length(T_noControl_pitch)
    Uout_noControl_pitch(i,:) = control_input0';
end
for i = 1:length(T_Control_pitch)
    de_nonlinear(i,1) = PitchAttitudeControl(theta_c, Y_Control_pitch(i,5), Y_Control_pitch(i,11), Kth, Kq);
    Uout_nonlinear_pitch(i,:) = control_input0' + [de_nonlinear(i) 0 0 0];
end


%% Part 3a) Roll Rate Control
%%% Transfer function from aileron to roll rate
% Using pure roll approximation phi_dot = delta p
% x_lat = [v, p, r, phi, psi, y] 


% [num_p2aileon, den_p2aileon] = ss2tf(Alat(1:4,1:4), Blat(1:4,1), [0 1 0 0],0);
[num_p2aileon, den_p2aileon] = ss2tf(Alat, Blat(:,1), [0 1 0 0 0 0],0);

num_roll_innerC = [0 -Ka 0 0 0 0];
den_roll_innerC = 1;

roll_innerC = feedback(tf(conv(num_roll_innerC,num_p2aileon), conv(den_roll_innerC,den_p2aileon)),1);
[Linear_num_roll_ss,Linear_den_roll_ss] = tfdata(roll_innerC,'v');
Nonlinear_roll_ss = ss(Alat,Blat(:,1),[0 1 0 0 0 0],0);
[wn_p_nonlinear,zeta_p_nonlinear] = damp(Nonlinear_roll_ss);

figure(102)
subplot(1,2,1) 
rlocus(num_p2aileon,den_p2aileon)
sgrid
title('Root Locus Roll Rate')

subplot(1,2,2)
step(roll_innerC,50); hold on
step(Nonlinear_roll_ss,50); hold on 
title('Roll Controller (Roll Rate to Aileon)')
legend(sprintf('Controlled Model (Ka=%.2f)', Ka), 'No-Control Model')

K_innerC = [0 -Ka 0 0 0 0];
New_Alat_inner = Alat - Blat(:,1)*K_innerC;
[New_Evec_inner,New_Eval_inner] = eig(New_Alat_inner);

% Non-Control and Control Simulations
Evec_RollRate = New_Evec_inner(:,6);
Evec_RollRate = real(Evec_RollRate/Evec_RollRate(2))*deg2rad(5); 
init_RollRate12_c = [0;Evec_RollRate(6);0;Evec_RollRate(4);0;Evec_RollRate(5);0;Evec_RollRate(1);0;Evec_RollRate(2);0;Evec_RollRate(3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Full sim in ode45
aircraft_state0 = trim_state + init_Roll12;
aircraft_state0_p_c= trim_state + init_RollRate12_c; 
control_input0 = trim_input; % [de, da, dr, dt];
TSPAN = [0 200];

[T_noControl_p,Y_noControl_p] = ode45(@(t,y) AircraftEOM(t,y,control_input0,wind_inertial,aircraft_parameters),TSPAN,aircraft_state0);
[T_Control_p,Y_Control_p] = ode45(@(t,y) AircraftEOMControl(t,y,control_input0,wind_inertial,aircraft_parameters,0,0,Ka,0,0),TSPAN,aircraft_state0_p_c);
% RollControl(ka, pc, dp) = da    pc = roll rate command, dp = delta roll
% rate

for i = 1:length(T_noControl_p)
    Uout_noControl_p(i,:) = control_input0';
end
for i = 1:length(T_Control_p)
    da_nonlinear_p(i,1) = RollControl(Ka,0,phi_c,Y_Control_p(i,4),Y_Control_p(i,10)); % ka, kp,phi_c,phi,p)
    Uout_nonlinear_roll_p(i,:) = control_input0' + [0 da_nonlinear_p(i) 0 0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Roll Angle & Roll Rate Control
% %%% Design the outer loop controller that feeds back roll angle to commanded roll rate
% x_lat = [v, p, r, phi, psi, y] 
[num_phi2p, den_phi2p] = ss2tf(Alat, Blat(:,1), [0 1 0 1 0 0], 0);
num_roll_outterC = [0 -Ka 0 -Kp*Ka 0 0];
den_roll_outterC = 1;

roll_outterC = feedback(tf(conv(num_roll_outterC,num_phi2p),conv(den_roll_outterC,den_phi2p)),1);
Nonlinear_phi2p_ss = ss(Alat,Blat(:,1),[0 1 0 1 0 0], 0);

figure(101)
subplot(1,2,1) 
rlocus(num_phi2p,den_phi2p)
sgrid
title('Root Locus Roll Angle to Roll Rate')

subplot(1,2,2)
step(roll_outterC,50); hold on
step(Nonlinear_phi2p_ss,50); hold on
title('Roll Controller (Roll Angle to Roll Rate)')
legend(sprintf('Linear Model (Kp=%.2f)', Kp), 'Nonlinear Model')

K_outterC = [0 -Ka 0 -Ka*Kp 0 0];
New_Alat_outter = Alat - Blat(:,1)*K_outterC;
[New_Evec_outter,New_Eval_outter] = eig(New_Alat_outter); 

% Non-Control and Control Simulations
Evec_Phi = New_Evec_outter(:,3);
Evec_Phi = real(Evec_Phi/Evec_Phi(4)*deg2rad(5)); 
init_Phi12_c = [0;Evec_Phi(6);0;Evec_Phi(4);0;Evec_Phi(5);0;Evec_Phi(1);0;Evec_Phi(2);0;Evec_Phi(3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Full sim in ode45
aircraft_state0 = trim_state + init_Roll12;
aircraft_state0_phi_c= trim_state + init_Phi12_c; 
control_input0 = trim_input; % [de, da, dr, dt];
TSPAN = [0 200];

[T_noControl_phi,Y_noControl_phi] = ode45(@(t,y) AircraftEOM(t,y,control_input0,wind_inertial,aircraft_parameters),TSPAN,aircraft_state0);
[T_Control_phi,Y_Control_phi] = ode45(@(t,y) AircraftEOMControl(t,y,control_input0,wind_inertial,aircraft_parameters,0,0,Ka,Kp,0),TSPAN,aircraft_state0_phi_c);
% RollControl(ka, pc, dp) = da    pc = roll rate command, dp = delta roll
% rate

for i = 1:length(T_noControl_phi)
    Uout_noControl_phi(i,:) = control_input0';
end
for i = 1:length(T_Control_phi)
    da_nonlinear_phi(i,1) = RollControl(Ka,Kp,phi_c,Y_Control_phi(i,4),Y_Control_phi(i,10)); % ka, kp,phi_c,phi,p)
    Uout_nonlinear_roll_phi(i,:) = control_input0' + [0 da_nonlinear_phi(i) 0 0];
end

%% Roll Control (Roll Rate + Roll Angle + Yaw Rate)
% %%% Design the outer loop controller that feeds back roll angle to commanded roll rate
% x_lat = [v, p, r, phi, psi, y] 
% [num_roll, den_roll] = ss2tf(Alat, Blat, [0 1 1 1 0 0], 0);
% num_roll_C = [0 Ka Kr Kp 0 0];
% den_roll_C = 1;
% 
% roll_C = feedback(tf(conv(num_roll_C,num_roll),conv(den_roll_C,den_roll)),1);
% Nonlinear_roll_ss = ss(Alat,Blat(:,1),[0 1 1 1 0 0], 0);
% 
% figure(102)
% step(roll_C,50); hold on
% step(Nonlinear_roll_ss,50); hold on
% title('Roll Controller (Roll Angle to Roll Rate)')
% legend('Linear Model', 'Nonlinear Model')

% Non-Control and Control Simulations
K_aileon = [0 -Ka 0 -Kp*Ka 0 0];
K_radder = [0 0 Kr 0 0 0];
New_Alat_roll = Alat - (Blat(:,1)*K_aileon + Blat(:,2)*K_radder);
[New_Evec_roll,New_Eval_roll] = eig(New_Alat_roll); 

Evec_Roll = New_Evec_roll(:,6);
Evec_Roll = real(Evec_Roll/Evec_Roll(4)*deg2rad(5)); 

init_Roll12_c = [0;Evec_Roll(6);0;Evec_Roll(4);0;Evec_Roll(5);0;Evec_Roll(1);0;Evec_Roll(2);0;Evec_Roll(3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Full sim in ode45
aircraft_state0 = trim_state + init_Roll12;
aircraft_state0_c= trim_state + init_Roll12_c; 
control_input0 = trim_input; % [de, da, dr, dt];
TSPAN = [0 200];

[T_noControl,Y_noControl] = ode45(@(t,y) AircraftEOM(t,y,control_input0,wind_inertial,aircraft_parameters),TSPAN,aircraft_state0);
[T_Control,Y_Control] = ode45(@(t,y) AircraftEOMControl(t,y,control_input0,wind_inertial,aircraft_parameters,0,0,Ka,Kp,Kr),TSPAN,aircraft_state0_c);
% RollControl(ka, pc, dp) = da    pc = roll rate command, dp = delta roll
% rate

for i = 1:length(T_noControl)
    Uout_noControl(i,:) = control_input0';
end
for i = 1:length(T_Control)
%     de_nonlinear(i,1) = PitchAttitudeControl(theta_c, Y_Control(i,5), Y_Control(i,11), Kth, Kq);
    da_nonlinear(i,1) = RollControl(Ka,Kp,phi_c,Y_Control(i,4),Y_Control(i,10)); % ka, kp,phi_c,phi,p)
    dr_nonlinear(i,1) = Kr*Y_Control(i,12);
    delta_input(i,:) = [0 da_nonlinear(i) dr_nonlinear(i) 0];
    Uout_nonlinear_roll(i,:) = control_input0' + delta_input(i);
end

%% Longitudinal + Lateral Dynamic Control 
%%% Full sim in ode45
aircraft_state0 = trim_state + init_Roll12 + init_Pitch12;
aircraft_state0_all_c= trim_state + init_Roll12_c + init_Pitch12_c; 
control_input0 = trim_input; % [de, da, dr, dt];
TSPAN = [0 200];

[T_noControl_all,Y_noControl_all] = ode45(@(t,y) AircraftEOM(t,y,control_input0,wind_inertial,aircraft_parameters),TSPAN,aircraft_state0);
[T_Control_all,Y_Control_all] = ode45(@(t,y) AircraftEOMControl(t,y,control_input0,wind_inertial,aircraft_parameters,Kth,Kq,Ka,Kp,Kr),TSPAN,aircraft_state0_all_c);
% RollControl(ka, pc, dp) = da    pc = roll rate command, dp = delta roll
% rate

for i = 1:length(T_noControl_all)
    Uout_noControl_all(i,:) = control_input0';
end
for i = 1:length(T_Control_all)
    de_nonlinear_all(i,1) = PitchAttitudeControl(theta_c, Y_Control_all(i,5), Y_Control_all(i,11), Kth, Kq);
    da_nonlinear_all(i,1) = RollControl(Ka,Kp,phi_c,Y_Control_all(i,4),Y_Control_all(i,10)); % ka, kp,phi_c,phi,p)
    dr_nonlinear_all(i,1) = Kr*Y_Control_all(i,12);
    delta_input_all(i,:) = [de_nonlinear_all(i) da_nonlinear(i) dr_nonlinear(i) 0];
    Uout_nonlinear_all(i,:) = control_input0' + delta_input_all(i);
end

%% Plots 
% Pitch Attitude Control (Longtitudinal Dynamic Control)
maintitle1 = 'Longitudinal Dynamic Control';
PlotAircraftSim(T_noControl_pitch, Y_noControl_pitch, Uout_noControl_pitch,[0;0;0],maintitle1,1,'r');
PlotAircraftSim(T_Control_pitch, Y_Control_pitch, Uout_nonlinear_pitch,[0;0;0],maintitle1,1,'b');

% Only Yaw Damper Control 
maintitle2 = 'Yaw Damper Control';
PlotAircraftSim(T_noControl_y, Y_noControl_y, Uout_noControl_y,[0;0;0],maintitle2,2,'r');
PlotAircraftSim(T_Control_y, Y_Control_y, Uout_nonlinear_roll_y,[0;0;0],maintitle2,2,'b');

% Only Roll Rate Control 
maintitle3 = 'Roll Rate Control';
PlotAircraftSim(T_noControl_p, Y_noControl_p, Uout_noControl_p,[0;0;0],maintitle3,3,'r');
PlotAircraftSim(T_Control_p, Y_Control_p, Uout_nonlinear_roll_p,[0;0;0],maintitle3,3,'b');

% Roll Angle & Roll Rate Control 
maintitle4 = 'Roll Rate + Roll Angle Control';
PlotAircraftSim(T_noControl_phi, Y_noControl_phi, Uout_noControl_phi,[0;0;0],maintitle4,4,'r');
PlotAircraftSim(T_Control_phi, Y_Control_phi, Uout_nonlinear_roll_phi,[0;0;0],maintitle4,4,'b');

% % Lateral Dynamic Control (Yaw Damper + Roll Controll)
maintitle5 = 'Laterial Dynamic Control';
PlotAircraftSim(T_noControl, Y_noControl, Uout_noControl,[0;0;0],maintitle5,5,'r');
PlotAircraftSim(T_Control, Y_Control, Uout_nonlinear_roll,[0;0;0],maintitle5,5,'b');

% Longtitudinal + Lateral Dynamic Control 
maintitle6 = 'Longitudinal + Laterial Dynamic Control';
PlotAircraftSim(T_noControl_all, Y_noControl_all, Uout_noControl_all,[0;0;0],maintitle6,6,'r');
PlotAircraftSim(T_Control_all, Y_Control_all, Uout_nonlinear_all,[0;0;0],maintitle6,6,'b');

%% Tables 
fprintf('Table 1: Eigenvalues (Longitudinal Dynamics)\n')
rowName1 = ["Original", "With Control"];
colName1 = ["1", "2", "3", "4", "5", "6"];
Eval_lon = [originalVal(1,1),originalVal(2,2),originalVal(3,3),originalVal(4,4),originalVal(5,5),originalVal(6,6)];
Eval_pitch_control = [New_eigenVal(1,1),New_eigenVal(2,2),New_eigenVal(3,3),New_eigenVal(4,4),New_eigenVal(5,5),New_eigenVal(6,6)];
table1 = array2table([Eval_lon; Eval_pitch_control],"RowNames", rowName1, "VariableNames", colName1);
disp(table1)

fprintf('Table 2: Eigenvalues (Yaw Damper)\n')
rowName1 = ["Original" ,"With Control"];
colName1 = ["1", "2", "3", "4", "5", "6"];
Eval_lat = [Eval_lat(1,1),Eval_lat(2,2),Eval_lat(3,3),Eval_lat(4,4),Eval_lat(5,5),Eval_lat(6,6)];
Eval_yaw = [New_Eval_yaw(1,1),New_Eval_yaw(2,2),New_Eval_yaw(3,3),New_Eval_yaw(4,4),New_Eval_yaw(5,5),New_Eval_yaw(6,6)];
table1 = array2table([Eval_lat; Eval_yaw],"RowNames", rowName1, "VariableNames", colName1);
disp(table1)

fprintf('Table 3: Eigenvalues (Roll Rate)\n')
rowName1 = ["Original" ,"With Control"];
colName1 = ["1", "2", "3", "4", "5", "6"];
Eval_p = [New_Eval_inner(1,1),New_Eval_inner(2,2),New_Eval_inner(3,3),New_Eval_inner(4,4),New_Eval_inner(5,5),New_Eval_inner(6,6)];
table1 = array2table([Eval_lat; Eval_p],"RowNames", rowName1, "VariableNames", colName1);
disp(table1)

fprintf('Table 4: Eigenvalues (Roll Rate & Roll Angle)\n')
rowName1 = ["Original" ,"With Control"];
colName1 = ["1", "2", "3", "4", "5", "6"];
Eval_outter = [New_Eval_outter(1,1),New_Eval_outter(2,2),New_Eval_outter(3,3),New_Eval_outter(4,4),New_Eval_outter(5,5),New_Eval_outter(6,6)];
table1 = array2table([Eval_lat; Eval_outter],"RowNames", rowName1, "VariableNames", colName1);
disp(table1)

fprintf('Table 5: Eigenvalues (Roll Control)\n')
rowName1 = ["Original" ,"With Control"];
colName1 = ["1", "2", "3", "4", "5", "6"];
Eval_roll = [New_Eval_roll(1,1),New_Eval_roll(2,2),New_Eval_roll(3,3),New_Eval_roll(4,4),New_Eval_roll(5,5),New_Eval_roll(6,6)];
table1 = array2table([Eval_lat; Eval_roll],"RowNames", rowName1, "VariableNames", colName1);
disp(table1)




%% ASEN 3111 - Computational Assignment 2 - Main
% Goal: The purpose of this assignment is to understand the airfoil theory
% using the superposition of elementary flows and to analyze the results
% with different flow parameters on streamlines, equipotential lines,
% and pressure contours. And also in order to validate the number of panels
% required, I find related error in velocity and pressure.

% Assumptions:
% 1. Thin symmetric airfoil
% 2. 2D incompressible, inviscid, & irrotational flow
% 3. Uniform flow + Vortex

% Author: Hideyuki Nakanishi
% SID: 110285397
% Date: Oct 09, 2022

clc;
clear all;
close all;
t0 = clock;
%% Constants
c = 2;              % (m) Chord length
alpha = 12*pi/180;  % (rad) Angle of attack = 12˚
V_inf = 68;         % (m/s) Freestream flow velocity
p_inf = 101.3*10^3; % (Pa) Freestream pressure
rho_inf = 1.225;    % (kg/m^3) Freestream density
N_true = 300;       % Number of discrete vortices for exact solution
N = 100;            % Number of discrete vortices
N_var = 1:N;        % Vector of different numbers of vorticies

% Chord Line
chord_x = linspace(0,c,N);
chord_y = zeros(1,N);


%% Exact Solutions (N_true Panels)
[StreamFunction, Phi, P_true, V_true, x, y] = Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N_true);
P_true_median = median(median(P_true)); % Median value of true pressure
V_true_median = median(median(V_true)); % Median value of true velocity 

%% Error Calculation
a = 0;  % Intialize it to count each blocks inside the P_error & V_error matrices
row = 0.5*length(N_var)*(1+length(N_var));  % This should be the length of the matrix
P_error = zeros(row,row);           % Initialize the P_error matrix
V_error = zeros(row,row);           % Initialize the V_error matrix
error_P = zeros(length(N_var),1);   % Initialize the error vector for pressure
error_V = zeros(length(N_var),1);   % Initialize the error vector for velocity

for i = N_var
    a = a + i;  % Count each blocks inside the matricies
    % Call the function to get pressure and velocity values corresponding
    % to N pannels 
    [~, ~, P_n, V_n, ~, ~] = Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,i);
    P_n(isnan(P_n)) = 0;            % Change NaN to 0
    V_n(isnan(V_n)) = 0;            % Change NaN to 0
    P_error(a-i+1:a,a-i+1:a) = P_n; % Save each P_n corresponding to N in the matrix
    V_error(a-i+1:a,a-i+1:a) = V_n; % Save each V_n corresponding to N in the matrix
    
    % Error calculation for pressure and velocity
    error_P(i) = abs(median(median(P_error(a-i+1:a,a-i+1:a))) - P_true_median)/P_true_median;
    error_V(i) = abs(median(median(V_error(a-i+1:a,a-i+1:a))) - V_true_median)/V_true_median;
    
end

%% Change in Each Variable (Chord Length, Angle of Attack, Freestream Velocity)
% Define the range of change in c, alpha, V_inf
n = 5;  % Number of the input values
change_c = [1,2,3,4,5];                   % (m) Change in chord length
change_a = [-20,-10,0,10,20].*pi./180;   % (rad) Change in angle of attack
change_v = [10,30,50,7,90];            % (m/s) Change in freestream velocity

for j = 1:n
    % Call the function to get outputs corresponding to each change in each variable
    [StreamFunction_c(N*(j-1)+1:N*j,:), Phi_c(N*(j-1)+1:N*j,:), ~, ~, x_c(N*(j-1)+1:N*j,:), y_c(N*(j-1)+1:N*j,:)] = Plot_Airfoil_Flow(change_c(j),alpha,V_inf,p_inf,rho_inf,N);
    [StreamFunction_a(N*(j-1)+1:N*j,:), Phi_a(N*(j-1)+1:N*j,:), ~, ~, x_a(N*(j-1)+1:N*j,:), y_a(N*(j-1)+1:N*j,:)] = Plot_Airfoil_Flow(c,change_a(j),V_inf,p_inf,rho_inf,N);
    [StreamFunction_v(N*(j-1)+1:N*j,:), Phi_v(N*(j-1)+1:N*j,:), ~, ~, x_v(N*(j-1)+1:N*j,:), y_v(N*(j-1)+1:N*j,:)] = Plot_Airfoil_Flow(c,alpha,change_v(j),p_inf,rho_inf,N);
    
    % Determine color levels for "contourf"
    levmin_psi_c = min(min(StreamFunction_c(N*(j-1)+1:N*j,:)));
    levmax_psi_c = max(max(StreamFunction_c(N*(j-1)+1:N*j,:)));
    levels_psi_c(j,:) = linspace(levmin_psi_c,levmax_psi_c,50);
    
    levmin_phi_c = min(min(Phi_c(N*(j-1)+1:N*j,:)));
    levmax_phi_c = max(max(Phi_c(N*(j-1)+1:N*j,:)));
    levels_phi_c(j,:) = linspace(levmin_phi_c,levmax_phi_c,50);
    
    levmin_psi_a = min(min(StreamFunction_a(N*(j-1)+1:N*j,:)));
    levmax_psi_a = max(max(StreamFunction_a(N*(j-1)+1:N*j,:)));
    levels_psi_a(j,:) = linspace(levmin_psi_a,levmax_psi_a,50);
    
    levmin_phi_a = min(min(Phi_a(N*(j-1)+1:N*j,:)));
    levmax_phi_a = max(max(Phi_a(N*(j-1)+1:N*j,:)));
    levels_phi_a(j,:) = linspace(levmin_phi_a,levmax_phi_a,50);
    
    levmin_psi_v = min(min(StreamFunction_v(N*(j-1)+1:N*j,:)));
    levmax_psi_v = max(max(StreamFunction_v(N*(j-1)+1:N*j,:)));
    levels_psi_v(j,:) = linspace(levmin_psi_v,levmax_psi_v,50);
    
    levmin_phi_v = min(min(Phi_v(N*(j-1)+1:N*j,:)));
    levmax_phi_v = max(max(Phi_v(N*(j-1)+1:N*j,:)));
    levels_phi_v(j,:) = linspace(levmin_phi_v,levmax_phi_v,50);
end

%% Plots
% Determine color levels for "contourf"
levmin_psi = min(min(StreamFunction));
levmax_psi = max(max(StreamFunction));
levels_psi = linspace(levmin_psi,levmax_psi,40);

levmin_phi = min(min(Phi));
levmax_phi = max(max(Phi));
levels_phi = linspace(levmin_phi,levmax_phi,40);

levmin_p = min(min(P_true));
levmax_p = max(max(P_true));
levels_p = linspace(levmin_p,levmax_p,40);

% Fig 1. Streamlines with N_true = 300 (True case)
figure(1)
hold on
contourf(x,y,StreamFunction,levels_psi); hold on
plot(chord_x(:).*cos(alpha),-chord_x(:).*sin(alpha),'LineWidth',2, 'color', 'r'); hold on
title(sprintf('Streamlines with N = %d', N_true))
xlabel('x')
ylabel('z')
colorbar
hold off

% Fig 2. Equipotential lines with N_true = 300 (True case)
figure(2)
hold on
contourf(x,y,Phi,levels_phi)
plot(chord_x.*cos(alpha),-chord_x.*sin(alpha),'LineWidth',2, 'color', 'r'); hold on
title(sprintf('Equipotential with N = %d', N_true))
xlabel('x')
ylabel('z')
colorbar
hold off

% Fig 3. Pressure Distribution with N_true = 300 (True case)
figure(3)
hold on
contourf(x,y,P_true,levels_p)
plot(chord_x.*cos(alpha),-chord_x.*sin(alpha),'LineWidth',2, 'color', 'r'); hold on
title(sprintf('Pressure Distribution with N = %d', N_true))
xlabel('x')
ylabel('z')
colorbar
hold off

% Error Associated with Pressure  & Velocity
figure(4)
hold on
subplot(1,2,1)
plot(N_var(2:end), error_P(2:end))
title('Pressure Field Accuracy Corresponding to N Panels')
ylabel('Error Associated with Pressure Field')
xlabel('Number of Panels')
grid on
hold off

hold on
subplot(1,2,2)
plot(N_var(2:end), error_V(2:end))
title('Resulting Flow Accuracy Corresponding to N Panels')
ylabel('Error Associated with Velocity')
xlabel('Number of Panels')
grid on
hold off

% Fig 5. Streamlines and Equipotential with Change in Chord Length
figure(5)
hold on
for k = 1:n
    subplot(2,n,k)
    contourf(x_c(N*(k-1)+1:N*k,:),y_c(N*(k-1)+1:N*k,:),StreamFunction_c(N*(k-1)+1:N*k,:),levels_psi_c(k,:))
    sgtitle('Streamlines (above) & Equipotential (below) with Change in Chord Length')
    title(sprintf('Chord  of %d m', change_c(k)))
    xlabel('x')
    ylabel('z')
    colorbar
    hold off
    hold on
    c_x = linspace(0,change_c(k),N);
    plot(c_x*cos(alpha),-c_x*sin(alpha),'LineWidth',2, 'color', 'r')
    hold off
    
    hold on
    subplot(2,n,k+n)
    contourf(x_c(N*(k-1)+1:N*k,:),y_c(N*(k-1)+1:N*k,:),Phi_c(N*(k-1)+1:N*k,:),levels_phi_c(k,:))
    title(sprintf('Chord of %d m', change_c(k)))
    xlabel('x')
    ylabel('z')
    colorbar
    hold off
    hold on
    plot(c_x*cos(alpha),-c_x*sin(alpha),'LineWidth',2, 'color', 'r')
    hold off
end

% Fig 6. Streamlines and Equipotential with Change in Angle of Attack
figure(6)
hold on
angle = [-20,-10,0,10,20];
for k = 1:n
    subplot(2,n,k)
    contourf(x_a(N*(k-1)+1:N*k,:),y_a(N*(k-1)+1:N*k,:),StreamFunction_a(N*(k-1)+1:N*k,:),levels_psi_a(k,:))
    sgtitle('Streamlines (above) & Equipotential (below) with Change in Angle of Attack')
    title(sprintf('Angle of %d˚', angle(k)))
    xlabel('x')
    ylabel('z')
    colorbar
    hold off
    hold on
    plot(chord_x*cos(change_a(k)),-chord_x*sin(change_a(k)),'LineWidth',2, 'color', 'r')
    hold off
    
    hold on
    subplot(2,n,k+n)
    contourf(x_a(N*(k-1)+1:N*k,:),y_a(N*(k-1)+1:N*k,:),Phi_a(N*(k-1)+1:N*k,:),levels_phi_a(k,:))
    title(sprintf('Angle of %d˚', angle(k)))
    xlabel('x')
    ylabel('z')
    colorbar
    hold off
    hold on
    plot(chord_x*cos(change_a(k)),-chord_x*sin(change_a(k)),'LineWidth',2, 'color', 'r')
    hold off
end

% Fig 7. Streamlines and Equipotential with Change in Velocity
figure(7)
hold on
for k = 1:n
    subplot(2,n,k)
    contourf(x_v(N*(k-1)+1:N*k,:),y_v(N*(k-1)+1:N*k,:),StreamFunction_v(N*(k-1)+1:N*k,:),levels_psi_v(k,:))
    sgtitle('Streamlines (above) & Equipotential (below) with Change in Freestream Velocity')
    title(sprintf('Velocity of %d m/s', change_v(k)))
    xlabel('x')
    ylabel('z')
    colorbar
    hold off
    hold on
    plot(chord_x*cos(alpha),-chord_x*sin(alpha),'LineWidth',2, 'color', 'r')
    hold off
    
    hold on
    subplot(2,n,k+n)
    contourf(x_v(N*(k-1)+1:N*k,:),y_v(N*(k-1)+1:N*k,:),Phi_v(N*(k-1)+1:N*k,:),levels_phi_v(k,:))
    title(sprintf('Velocity of %d m/s', change_v(k)))
    xlabel('x')
    ylabel('z')
    colorbar
    hold off
    hold on
    plot(chord_x*cos(alpha),-chord_x*sin(alpha),'LineWidth',2, 'color', 'r')
    hold off
end

time = etime(clock,t0);
fprintf('Time = %.1f s\n', time);



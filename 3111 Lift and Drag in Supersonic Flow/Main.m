%% ASEN 3111 - Computational Assignment 4 - Main
%
% Goal: Understanding the Shock-Expanssion Theory and Linearlized Theory.
% In addition, To understand the impact of change in Mach and angle of
% attack on lift and wave drag coefficients
%
% Author: Hideyuki Nakanishi
% SID: 110285397
% Date: Dec 04, 2022

clc
clear
close all
t0 = clock;

%% Q1 Validation of the ø - ß - M Diagram
Gamma = 1.4;            % Gamma in the thermally calorically perfect condition
N = 1000;               % Number of theta
Thetad = linspace(1,50,N);  % Vector of wedge angles
% Chosing Mach numbers
M1 = 1.1:0.05:1.5;
M2 = 1.6:0.1:2.0;
M3 = 2.2:.2:4.0;
M4 = 4.5;
M5 = 5:1:10;
M_vec = [M1, M2, M3, M4, M5];   % Vector of Mach numbers

% Call the 'ObliqueShockBeta' function for weak and strong to get beta
for i = 1:length(M_vec)
    count_w = 0;    % Count the real postive weak solutions
    count_s = 0;    % Count the real postive strong solutions
    for j = 1:length(Thetad)
        Beta_weak(i,j) = ObliqueShockBeta(M_vec(i),Thetad(j),Gamma,'Weak');     % Find the weak solutions
        Beta_strong(i,j) = ObliqueShockBeta(M_vec(i),Thetad(j),Gamma,'Strong'); % Find the strong solutions
        
        if Beta_weak(i,j) >= 0 && isreal(Beta_weak(i,j))  % Ignore non-real, negative numbers
            count_w = count_w + 1;
            Thetad_w(i,count_w) = Thetad(j);
            Beta_w(i,count_w) = Beta_weak(i,j);
        end
        if Beta_strong(i,j) >= 0 && isreal(Beta_strong(i,j))
            count_s = count_s + 1;
            Thetad_s(i,count_s) = Thetad(j);
            Beta_s(i,count_s) = Beta_strong(i,j);
        end
    end
    % Beta matrix (Row: Mach number Col: Theta)
    Beta = [Beta_w(i,:), fliplr(Beta_s(i,1:count_s-1))];
    Thetad = [Thetad_w(i,:), fliplr(Thetad_s(i,1:count_s-1))];
    % Plot the theta-beta-m relations
    figure(1)
    hold on; axis tight; grid on
    plot(Thetad,Beta);
    title('θ − β − M Diagram');
    ylabel('Beta β (deg)');
    xlabel('Theta θ (deg)');
    Thetad = linspace(1,50,N);  % Reset it to the original thetad vector
end


%% Q2 Computation of Lift and Drag for a Diamond-Wedge Airfoil
M = 2;          % Freestream Mach number
alpha = 5;      % [deg] Angle of attack
epsilon1 = 5;  % [deg] Geometric angle of diamond airfoil in front
epsilon2 = 10;   % [deg] Geometric angle of diamond airfoil in back

% Get sectional lift and wave drag coefficients by calling 'DiamondAirfoil'
[c_l, cd_w] = DiamondAirfoil(M,alpha,epsilon1, epsilon2);
fprintf('Lift Coefficient = %.3f\n', c_l);
fprintf('Wave Drag Coefficient = %.3f\n', cd_w);

%% Q3 Impact of Angle of Attack and Mach Number
n = 200;
alpha = linspace(-25,25,n); % [deg] Angle of attack
M = [2,3,4,5];  % Mach numbers
c = 3;          % [m] Chord length

% Initialization
cl_linear = zeros(n,length(M));   % Linearized lift coefficient
cd_w_linear = zeros(n,length(M)); % Linearized wave drag coefficient
cl = zeros(n,length(M));          % Lift coefficient from Shock-Expanssion 
cd_w = zeros(n,length(M));        % Wave drag coefficient from Shock-Expanssion

% Get lift and drag coeff with different Mach over angle of attack
for i = 1:length(M) 
    [cl_linear(:,i),cd_w_linear(:,i)] = LinearizedTheory(c,alpha,M(i),epsilon1,epsilon2);
    [cl(:,i),cd_w(:,i)] = DiamondAirfoil(M(i),alpha,epsilon1,epsilon2);
    
    % Plots: Lift Coeff: Linearized Theory vs. Shock-Expanssion Theory
    figure(2)
    clr = ['r','b','g','k'];
    p1(i) = plot(alpha, cl_linear(:,i), 'DisplayName', sprintf('Linearized Theory (M = %d)', M(i))); hold on
    p2(i) = plot(alpha, cl(:,i), '--','DisplayName', sprintf('Shock-Expanssion (M = %d)',M(i)));
    title('Cl - Linearized Theory vs. Shock-Expanssion Theory over Angle of Attack');
    xlabel('Angle of Attack (deg)')
    ylabel('Lift Coefficient')
    set([p1(i),p2(i)],'Color',clr(i))
    grid on
    
    % Plots: Wave Drag Coeff: Linearized Theory vs. Shock-Expanssion Theory
    figure(3)
    p3(i) = plot(alpha, cd_w_linear(:,i), 'DisplayName', sprintf('Linearized Theory (M = %d)', M(i))); hold on
    p4(i) = plot(alpha, cd_w(:,i), '--','DisplayName', sprintf('Shock-Expanssion (M = %d)',M(i)));
    title('Cd - Linearized Theory vs. Shock-Expanssion Theory over Angle of Attack');
    xlabel('Angle of Attack (deg)')
    ylabel('Wave Drag Coefficient')
    set([p3(i),p4(i)],'Color',clr(i))
    grid on
end
legend([p1,p2],'Location', 'north')
legend([p3,p4],'Location', 'north')
legend show


%% Time Taken
time = etime(clock,t0);
fprintf('Time = %.1f s\n', time);
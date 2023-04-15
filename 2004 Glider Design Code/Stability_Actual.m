% Assumptions *********************************
% - SLUF ( L = W )
% - STD ATM @ 1500m

% Requirements ********************************
R_req = [125,175];      % (m) Min - Max of glider range
V_Rreq = [7, 12];       % (m/s) minimum velocity for range
E_req = [13, 20];       % (s) maximum time of flight
theta_req = [10, 20];   % (deg) Elevetor pitch control
vol_HT_req = [0.3,0.6];     % Volume coefficient of horizontal tail
vol_VT_req = [0.02,0.05];   % Volume coefficient of vertical tail
SM_req = [0.15,0.25];       % Longitudinal Stability

% Constants ***********************************
chord = 0.18;       % Chord length of wing body
b = 1;              % (m) Length of wing span
S_ref = b*chord;    % Surface area of wing body
x_h = 0.101;       
S_ht = 0.036;
b_h = S_ht/x_h;
x_v = x_h;
S_vt = 2*0.007;
S_f = 0.033;
S_wet = 2*S_ref + 2*S_ht + 2*S_vt + S_f;        % Wetted area of aircraft

e0 = 0.7;           % Efficiency factor
h = 17.5;           % (m) Height where a glider is thrown at
rho = 1.0581;       % (kg/m^3) At an altitude of 1500m
C_fe = 0.003;       % Equivalent coefficient of skin friction
ws = 0.295*9.81;    % (N/m^2) Dollar Tree Foam's approximate specific weight
W_pay = 0.16*9.81;  % (N) Weigth of camera

d_tail_dwn = 0.4;   % Coefficient of the impact of tail downwash
CM_ac = 0;       % Coefficient of moment relative to the aerodynamic center of a flat plate
B = 5;              % Spiral parameter: Neutral B = 5; Stable B > 5

i_t = 0;            % (deg) Angle of incident 
[a,CL_w] = slp.calc_slope(S_wet,S_ref,b);       % Lift slope and CL for wing
[a_t,CL_tail] = slp.calc_slope(S_wet,S_ht,b_h); % Lift slope and CL fpr tail
e0_dwn = 0;         % Downwash effect at zero-lift
alpha_w = CL_w/a;   % (deg) Angle of attack for wing body
e_dwn = e0_dwn + d_tail_dwn*alpha_w;            % Downwash effect
alpha_t = alpha_w - i_t - e_dwn;                % (deg) Angle of attack for tail

x0 = 0.137;     % (m) Distance FROM the tip of NOSE TO the tip of WING body
l_cg = 0.301;   % (m) Position of cg FROM the tip of the NOSE
l_f = 0.68;     % (m) Length of the fuselage FROM the tip of the NOSE
x_cg = l_cg - x0;   % (m) Position of cg FROM the tip of the WING
x_ac_w = chord/4;   % Aerodynamic center of wing body: quater chord

% Horizontal Tail 
ac_h = x_h/4;
x_ac_h = l_f - x0 - (x_h - ac_h);
vol_H = S_ht*(x_ac_h - x_cg)/(S_ref*chord);
CL_t = a_t*alpha_w*(1-d_tail_dwn) - a_t*(i_t + e0_dwn);

% Vertical Tail 
x_ac_v = x_ac_h;
vol_V = S_vt*(x_ac_v - x_cg)/(S_ref*b);

% Moment Coefficient
h_ac_w = x_ac_w/chord;
h_cg = x_cg/chord;
h_NP = h_ac_w + vol_H*(a_t/a)*(1 - d_tail_dwn);
SM = h_NP - h_cg; % Should be positive 
d_CM_cg = -a*SM;
CM0 = CM_ac + vol_H*a_t*(i_t + e0_dwn);
CM_cg = CM0 + alpha_w*d_CM_cg;

fprintf('          Tail Coefficient Lift: %.1f\n', CL_t);
fprintf('           Tail Angle of Attack: %.1f\n', alpha_t);
fprintf('  Horizontal Tail Volume (V_H) : %.1f\n', vol_H);
fprintf('    Vertical Tail Volume (V_v) : %.3f\n', vol_V);
fprintf('  Coefficient of Lift for Tail : %.3f\n', CL_t);
fprintf('                CM at zerolift : %.3f\n', CM0);
fprintf('CM around the Center of Gravity: %.3f\n',CM_cg);
fprintf('            Static Margin (SM) : %.1f percent\n', abs(SM)*100);

figure(3)
FP_Langle = [-13.0728 ,-12.09 ,-11.0735 ,-10.0574,-9.07559 ,-8.09424 ,-7.14725 ,-6.13293,-5.05147 ,-4.07163 ,-3.09146 ,-2.07747,0.930398 ,1.97819 ,2.99231 ,4.00636 ,5.05408 ,7.0479 ,8.02842 ,9.04358 ,10.0595 ,11.0079 ,12.0582,13.0068];
FP_Cl = [-0.787245,-0.804155,-0.815366,-0.806654,-0.783716,-0.740854,-0.678074,-0.598204,-0.492707,-0.387227,-0.295978,-0.201877,0.0889578,0.185911,0.274319,0.365573,0.465372,0.664953,0.74197,0.787684,0.802089,0.805097,0.799584,0.794054];
Cl_SlopeLine = polyfit(FP_Langle(5:19),FP_Cl(5:19),1); %polyfit output slope and y intercept
Cl_Slope = @(x) Cl_SlopeLine(1)*x + Cl_SlopeLine(2);%y= mx+b           
alpha0 = (-1)*Cl_SlopeLine(2)/Cl_SlopeLine(1);%finding alpha at 0
wing_3dCl = a * (FP_Langle(:) - alpha0);%final 3d wing Cl calculation
tail_3dCl = a_t * (FP_Langle(:) - alpha0);
plot(FP_Langle, FP_Cl)
hold on
plot(FP_Langle, Cl_Slope(FP_Langle))%best fit
hold on
plot(FP_Langle, wing_3dCl);
plot(FP_Langle, tail_3dCl);
title('Angle of Attack vs Lift Coefficient');
xlabel('Angle of Attack');
ylabel('Coefficient of Lift');
legend('2d truth', 'best fit', 'Wing Lift Curve', 'Horizontal Tail Lift Curve');
grid on
hold off

figure(6)
hold on 
FP_Langle = [-13.0728 ,-12.09 ,-11.0735 ,-10.0574,-9.07559 ,-8.09424 ,-7.14725 ,-6.13293,-5.05147 ,-4.07163 ,-3.09146 ,-2.07747,0.930398 ,1.97819 ,2.99231 ,4.00636 ,5.05408 ,7.0479 ,8.02842 ,9.04358 ,10.0595 ,11.0079 ,12.0582,13.0068];
y = @(x) d_CM_cg*x + CM0;
trim = -CM0/d_CM_cg;
plot(FP_Langle,y(FP_Langle));
xline(trim,'--r',{'Trimmed'});
title('Trim Diagram Based on Our Actual Design');
ylabel('∂CMcg/∂α');
xlabel('α_w');
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
xlim([-2.5 10])
ylim([-0.2 0.2])
hold off


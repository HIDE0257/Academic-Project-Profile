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
b_h = 0.35;         % Wingspan of horizontal tail
S_ref = b*chord;    % Reference area of wing body
h = 17.5;           % (m) Height where a glider is thrown at
C_fe = 0.003;       % Equivalent coefficient of skin friction
rho = 1.0581;       % (kg/m^3) At an altitude of 1500m
e0 = 0.7;           % Efficiency factor
B = 5;              % Spiral parameter: Neutral B = 5; Stable B > 5
vol_v = 0.035;      % Desired vertical tail vol coefficient
ws = 0.295*9.81;    % (N/m^2) Dollar Tree Foam's approximate specific weight
W_pay = 0.16*9.81;  % (N) Weigth of camera
CM_ac = 0;       % Coefficient of moment relative to the aerodynamic center of a flat plate
d_tail_dwn = 0.4;   % Coefficient of the impact of tail downwash
x_ac_w = chord/4;   % Aerodynamic center of wing body: quater chord
e0_dwn = 0;         % Downwash effect at zero-lift

x0 = 0.137;     % (m) Distance from the tip of nose to the tip of wing body
x_cg = 0.03926;  % (m) Position of cg from the tip of the wing
S_f = 0.033;

S_ht = linspace(0.01,0.5,100);   % (m^2) Horizontal tail surface area
S_vt = linspace(0.001,0.5,100);
x_ac_h = linspace(0.3,1,100);
i_t = linspace(0.1,20,100);       % (deg) Incident angle of tail

h_ac_w = x_ac_w/chord;
h_cg = x_cg/chord;

for i = 1:length(S_ht)
    for j = 1:length(S_vt)
        for k = 1:length(x_ac_h)
            % Center of Gravity 
            x_h(i) = S_ht(i)/b_h;
            x_v(j) = x_h(i);
            % Horizontal Tail
            ac_h(i) = x_h(i)/4;
            vol_H(i,k) = S_ht(i)*(x_ac_h(k) - x_cg)/(S_ref*chord);
            
            
            % Vertical Tail
            vol_V(j,k) = S_vt(j)*(x_ac_h(k) - x_cg)/(S_ref*b);
            
            if vol_HT_req(1) <= vol_H(i,k) && vol_H(i,k) <= vol_HT_req(2)
                dim_ht(:,i) = [x_h(i), b_h, S_ht(i)];
                VH(i,k) = vol_H(i,k);
                L_f(i,k) = x_ac_h(k) + x0 + 3*x_h(i)/4;
                x_ac_h(k) = x_ac_h(k);
            end
            
            if vol_VT_req(1) <= vol_V(j,k) && vol_V(j,k) <= vol_VT_req(2)
                dim_vt(j) = S_vt(j);
                VV(j,k) = vol_V(j,k);
            end
        end
    end
end

% Pick up the values from the array
[~,c_VH] = find(dim_ht(1,:) > 0);
S_vt = dim_vt(1,3);
vol_V = VV(3,1);
n = 1;
for p = c_VH
    for q = 1:length(i_t)
        Swet(p) = 2*S_ref + 2*dim_ht(3,p) + 2*S_vt + S_f;
        [a(p),CL_w(p)] = slp.calc_slope(Swet(p),S_ref,b);
        [a_t(p),cl_t(p)] = slp.calc_slope(Swet(p),dim_ht(3,p),b_h);
        alpha_w(p) = CL_w(p)/a(p);
        e_dwn(p) = e0_dwn + d_tail_dwn*alpha_w(p);   % Downwash effect
        h_ac_w = x_ac_w/chord;
        h_NP(p) = h_ac_w + VH(p,1)*(a_t(p)/a(p))*(1 - d_tail_dwn);
        h_cg = x_cg/chord;
        SM(p) = h_NP(p) - h_cg;
        d_CM_cg(p) = -a(p)*SM(p);
        
        alpha_t(p,q) = alpha_w(p) - i_t(q) - e_dwn(p);
        CL_t(p,q) = a_t(p)*alpha_w(p)*(1-d_tail_dwn) - a_t(p)*(i_t(q) + e0_dwn);
        % Moment Coefficient
        CM0(p,q) = CM_ac + VH(p,1)*a_t(p)*(i_t(q) + e0_dwn);
        CM_cg(p,q) = CM0(p,q) + alpha_w(p)*d_CM_cg(p);
        
        if CM_cg(p,q) < 0
            if CM0(p,q) > 0
                if SM_req(1) <= SM(p) && SM(p) <= SM_req(2)
                    finalCM_cg(p,q) = CM_cg(p,q);
                    finalCM0(p,q) = CM0(p,q);
                    finali_t(q) = i_t(q);
                    finalSwet(p) = Swet(p);
                    finalCL_t(p,q) = CL_t(p,q);
                    finalS_h(p) = dim_ht(3,p);
                    finalL_F(p,1) = L_f(p,1);
                    finalSM(p) = SM(p);
                    finalx_ac_h(p) = x_ac_h(p);
                    a(p) = a(p);
                end
            end
        end
    end
end

row = find(SM > 0);
col = find(i_t > 0);
c = row(7);
r = col(14);

trim = -CM0(row,col)/d_CM_cg(c);
fprintf('                  Wetted Area: %.3f m^2\n', finalSwet(c));
fprintf('                       x_ac,h: %.2f m\n', finalx_ac_h(c));
fprintf('           Length of Fuselage: %.2f m\n', finalL_F(c));
fprintf('     Hori Tail Reference area: %.3f m^2\n', finalS_h(c));
fprintf('     Vert Tail Reference Area: %.3f m^2\n', S_vt);
fprintf('Coefficient of Lift for Tail : %.3f\n', finalCL_t(r,c));
fprintf('                         Trim: %.1f\n', trim);
fprintf('               CM at zerolift: %.3f\n', CM0(r,c));
fprintf('CM around the Center of Gravity: %.3f\n', CM_cg(r,c));
fprintf('               Incident Angle: %.2f\n', finali_t(c));
fprintf('                Static Margin: %.1f percent\n', finalSM(c)*100);

S_wet = 2*S_ref + 2*finalS_h(c) + 2*S_vt + S_f;
W_TO = ws*S_wet + W_pay;
M_w = (chord/2)*(ws*S_ref);
M_f = ((finalL_F(b)/2)-x0)*(ws*S_f);
M_ht = (finalx_ac_h(b)+(ac_h(c)/4))*(ws*finalS_h(c));
M_vt = (finalx_ac_h(b)+(ac_h(c)/4))*(ws*S_vt);
x_cg = (M_w+M_f+M_ht+M_vt)/W_TO;
fprintf(' Centre of Gravity: %.5f m\n', x_cg);

figure(5)
FP_Langle = [-13.0728 ,-12.09 ,-11.0735 ,-10.0574,-9.07559 ,-8.09424 ,-7.14725 ,-6.13293,-5.05147 ,-4.07163 ,-3.09146 ,-2.07747,0.930398 ,1.97819 ,2.99231 ,4.00636 ,5.05408 ,7.0479 ,8.02842 ,9.04358 ,10.0595 ,11.0079 ,12.0582,13.0068];
slope = @(x) d_CM_cg(c)*x + CM0(row,col);
hold on
plot(FP_Langle,slope(FP_Langle));
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

% rowName = ["Tail Chord Legth (m)", "Tail Surface Area (m^2)", "Vertical Tail Height", "Vertical Tail Area (m^2)", "Pitching Moment", "Angle", "Incident Angle", "Dihedral Angle"];
% best_val = array2table(BEST,"RowNames", rowName); % "VariableNames", varName);

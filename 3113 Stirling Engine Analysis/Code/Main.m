% Thermodynamics Lab - 02/28/2023
clc;
clear all;
close all;

%% Constants
voltage = 48;   %[V] voltage for the heater
R = 287;        %[J/kg*K] Ideal Gas constant
rho = 1.14;     %[m^3] Air density at an altitude of Boulder
V_cylinder = 1.7267e-4; %[kg/m^3] Volume of the cylinder

%% Extracting Data 
%8 deg
data_8 = load("8degGroupEight");
time_8 = data_8(:,1);
pressure_8 = data_8(:,2);
temp_topoftop_8 = data_8(:,3);
temp_bottomoftop_8 = data_8(:,4);
temp_topofbottom_8 = data_8(:,5);
temp_bottomofbottom_8 = data_8(:,6);
optical_switch_8 = data_8(:,8);

%10 deg
data_10 = load("10degGroupEight");
time_10 = data_10(:,1);
pressure_10 = data_10(:,2);
temp_topoftop_10 = data_10(:,3);
temp_bottomoftop_10 = data_10(:,4);
temp_topofbottom_10 = data_10(:,5);
temp_bottomofbottom_10 = data_10(:,6);
optical_switch_10 = data_10(:,8);

%12 deg
data_12 = load("12degGroupEight");
time_12 = data_12(:,1);
pressure_12 = data_12(:,2);
temp_topoftop_12 = data_12(:,3);
temp_bottomoftop_12 = data_12(:,4);
temp_topofbottom_12 = data_12(:,5);
temp_bottomofbottom_12 = data_12(:,6);
optical_switch_12 = data_12(:,8);

%% Finding RPM
%8 deg
points_8 = []; 
for i = 1:length(optical_switch_8) %look through whole data set
    if optical_switch_8(i) == 1 %check if optical switch is at true
        if last == 0 %check if last value of optical switch was true
            points_8 = [points_8; i]; %if optical switch was last false, record point
        end
        last = 1; %set last to true, which will save the last state for next iteration
    elseif optical_switch_8(i) == 0 %if optical switch is false
    last = 0; %set last to 0, so that the next point will be recorded
    end
end

time_steps_8 = time_8(points_8); %use the points collected to find the corresponding times
timespan_8 = (time_steps_8(end)-time_steps_8(1))/60; %find the time between first and last points, convert to minutes
avg_rpm_8 = (length(time_steps_8)-1)/timespan_8; %divide number of cycles by time span to find average rpm
rpm_8 = [];
for i = 1:length(time_steps_8)-1
    rpm_8 = [rpm_8; 1/((time_steps_8(i+1)-time_steps_8(i))/60)];
end

%10 deg
points_10 = [];
for i = 1:length(optical_switch_10)
    if optical_switch_10(i) == 1
        if last == 0
            points_10 = [points_10; i];
        end
        last = 1;
    elseif optical_switch_10(i) == 0
    last = 0;
    end
end

time_steps_10 = time_10(points_10);
timespan_10 = (time_steps_10(end)-time_steps_10(1))/60;
avg_rpm_10 = (length(time_steps_10)-1)/timespan_10;
rpm_10 = [];

for i = 1:length(time_steps_10)-1
    rpm_10 = [rpm_10; 1/((time_steps_10(i+1)-time_steps_10(i))/60)];
end

%12 deg
points_12 = []; 
for i = 1:length(optical_switch_12) %look through whole data set
    if optical_switch_12(i) == 1 %check if optical switch is at true
        if last == 0 %check if last value of optical switch was true
            points_12 = [points_12; i]; %if optical switch was last false, record point
        end
        last = 1; %set last to true, which will save the last state for next iteration
    elseif optical_switch_12(i) == 0 %if optical switch is false
    last = 0; %set last to 0, so that the next point will be recorded
    end
end

time_steps_12 = time_12(points_12); %use the points collected to find the corresponding times
timespan_12 = (time_steps_12(end)-time_steps_12(1))/60; %find the time between first and last points, convert to minutes
avg_rpm_12 = (length(time_steps_12)-1)/timespan_12; %divide number of cycles by time span to find average rpm

rpm_12 = [];
for i = 1:length(time_steps_12)-1
    rpm_12 = [rpm_12; 1/((time_steps_12(i+1)-time_steps_12(i))/60)];
end

%% PV Diagram (SolidWork Part) 
fgdjkl = 342006.34;
foam_volume = 169331.84;

volume_container = fgdjkl - foam_volume;
pp = 176.71; % [mm^2] Area of the piston

% [mm] Displacement of the piston (Stroke)
displacement_data1 = readmatrix("smll_btm_fce_disp_70.6.csv");
displacement_data2 = readmatrix("Small_Bottom_Face_Disp_98.4.csv");
displacement_data3 = readmatrix("Small_Bottom_Face_Disp2_116.csv");

% [mm^3] Change in volume by work 
volume_run1 = displacement_data1(:,2) .* pp;
volume_run2 = displacement_data2(:,2) .* pp;
volume_run3 = displacement_data3(:,2) .* pp;

% [m^3] Total Volume
v_total_1 = (volume_run1 + volume_container)*10^-9;
v_total_2 = (volume_run2 + volume_container)*10^-9;
v_total_3 = (volume_run3 + volume_container)*10^-9;

% [psi] Pressure data from the sensor
pressure_data1 = readmatrix("8degGroupEight");
pressure_data2 = readmatrix("10degGroupEight");
pressure_data3 = readmatrix("12degGroupEight");

% [Pa] Convert psi to Pa
psi1 = pressure_data1(:,2) .* 6895;
psi2 = pressure_data2(:,2) .* 6895;
psi3 = pressure_data3(:,2) .* 6895;

% [s] Time vector from the preesure data 
xq1 =  pressure_data1(:,1); 
xq2 =  pressure_data2(:,1);
xq3 =  pressure_data3(:,1);

% Interporation 
vq1 = interp1(xq1,psi1,displacement_data1(:,1));
vq2 = interp1(xq2,psi2,displacement_data2(:,1));
vq3 = interp1(xq3,psi3,displacement_data3(:,1));

% [Pressure & Volume] One complete cycle for each dataset
vq1 = -vq1(104:232);            % [Pa] Pressure 
v_total_1=v_total_1(104:232);   % [m^3] Volume
Pmax1 = max(vq1);               % [Pa] Maximum Pressure
    
vq2 = -vq2(96:218);             % [Pa] Pressure 
v_total_2=v_total_2(96:218);    % [m^3] Volume
Pmax2 = max(vq2);               % [Pa] Maximum Pressure

vq3 = -vq3(107:233);            % [Pa] Pressure  
v_total_3=v_total_3(107:233);   % [m^3] Volume
Pmax3 = max(vq3);               % [Pa] Maximum Pressure

%% PV Diagrams
% Data1 - 8 deg
figure
% subplot(3,1,1)
hold on
plot(v_total_1,vq1)
title('Data 1: PV-Diagram (8 deg)')
xlabel('Volume (m^3)')
ylabel('Pressure (Pa)')
grid
hold off

figure
% Data2 - 10 deg
% subplot(3,1,2)
hold on 
plot(v_total_2,vq2)
title('Data 2: PV-Diagram (10 deg)')
xlabel('Volume (m^3)')
ylabel('Pressure (Pa)')
grid
hold off

figure
% Data3 - 12 deg
% subplot(3,1,3)
hold on 
plot(v_total_3,vq3)
title('Data 3: PV-Diagram (12 deg)')
xlabel('Volume (m^3)')
ylabel('Pressure (Pa)')
grid
hold off


%% Current Data & Ideal Qin Calculation
A8 = load('current8deg.mat');
A10 = load('current10deg.mat');
A12 = load('current12deg.mat');

% % To know the average cycle time for each dataset 
% for i = 2:length(time_steps_8) 
%     dt_8(i-1) = time_steps_8(i) - time_steps_8(i-1);
% end
% for i = 2:length(time_steps_10) 
%     dt_10(i-1) = time_steps_10(i) - time_steps_10(i-1);
% end
% for i = 2:length(time_steps_8) 
%     dt_12(i-1) = time_steps_12(i) - time_steps_12(i-1);
% end
% avg_dt_8 = mean(dt_8);
% avg_dt_10 = mean(dt_10);
% avg_dt_12 = mean(dt_12);

% Area under the curve of current with respect to time
Idt_8 = trapz(A8.dat8(:,1),A8.dat8(:,2));
Idt_10 = trapz(A10.dat10(:,1),A10.dat10(:,2));
Idt_12 = trapz(A12.dat12(:,1),A12.dat12(:,2));

% Assume Qin = Win --> Electricity power (Win = V*Idt (J))
Qin_ideal_8 = voltage*Idt_8; 
Qin_ideal_10 = voltage*Idt_10; 
Qin_ideal_12 = voltage*Idt_12; 


%% Efficiency
%% Carnot Efficiency (n = 1 - T_L/T_H) (%)
T_H_8 = mean(temp_topofbottom_8(points_8(3):points_8(4))) + 273.15; %Finding average temperature for one cycle on bottom (T_low)
T_L_8 = mean(temp_bottomoftop_8(points_8(3):points_8(4))) + 273.15; %Finding average temperature for one cycle on top (T_high)
carnot_8 = (1 - T_L_8/T_H_8) * 100; %Calculating percentage efficiency

T_H_10 = mean(temp_topofbottom_10(points_10(3):points_10(4))) + 273.15; %Finding average temperature for one cycle on bottom (T_low)
T_L_10 = mean(temp_bottomoftop_10(points_10(3):points_10(4))) + 273.15; %Finding average temperature for one cycle on top (T_high)
carnot_10 = (1 - T_L_10/T_H_10) * 100; %Calculating percentage efficiency

T_H_12 = mean(temp_topofbottom_12(points_12(3):points_12(4))) + 273.15; %Finding average temperature for one cycle on bottom (T_low)
T_L_12 = mean(temp_bottomoftop_12(points_12(3):points_12(4))) + 273.15; %Finding average temperature for one cycle on top (T_high)
carnot_12 = (1 - T_L_12/T_H_12) * 100; %Calculating percentage efficiency

carnot = [carnot_8, carnot_10, carnot_12]';

%% Experimental Efficiency (%)
% Get the max and min volumes for each dataset
[Vmax_8,idx_max_8] = max(v_total_1);
[Vmin_8,idx_min_8] = min(v_total_1);

[Vmax_10,idx_max_10] = max(v_total_2);
[Vmin_10,idx_min_10] = min(v_total_2);

[Vmax_12,idx_max_12] = max(v_total_3);
[Vmin_12,idx_min_12] = min(v_total_3);

% Find two x-intersect points to calculate the area of the PV diagram
n = 0;
left_1 = [];
right_1 = [];
for i = 1:length(vq1)
    if sign(vq1(i)) == 1
        n = n + 1;
        left_1(n) = i;
    end
end
left_1 = left_1(1);
right_1 = left_1 + n;

n = 0;
left_2 = [];
right_2 = [];
for i = 1:length(vq2)
    if sign(vq2(i)) == 1
        n = n + 1;
        left_2(n) = i;
    end
end
left_2 = left_2(1);
right_2 = left_2 + n;

n = 0;
left_3 = [];
right_3 = [];
for i = 1:length(vq3)
    if sign(vq3(i)) == 1
        n = n + 1;
        left_3(n) = i;
    end
end
left_3 = left_3(1);
right_3 = left_3 + n;

% Get the areas of the top and bottom of the circle --> Wnet_out (J)
top_area_8 = abs(trapz(v_total_1(idx_min_8:right_1),vq1(idx_min_8:right_1))) - abs(trapz(v_total_1(left_1:idx_min_8),vq1(left_1:idx_min_8)));
bottom_area_8 = abs(trapz(v_total_1(idx_max_8:end),vq1(idx_max_8:end))) + abs(trapz(v_total_1(1:left_1),vq1(1:left_1))) - abs(trapz(v_total_1(right_1:idx_max_8),vq1(right_1:idx_max_8)));
Wnet_8 = top_area_8 + bottom_area_8;    % [J] Wnet_out for 8 deg
exp_8 = (Wnet_8/Qin_ideal_8)*100;             % [%] Actual efficiency for 8 deg

top_area_10 = abs(trapz(v_total_2(idx_min_10:right_2),vq2(idx_min_10:right_2))) - abs(trapz(v_total_2(left_2:idx_min_10),vq2(left_2:idx_min_10)));
bottom_area_10 = abs(trapz(v_total_2(idx_max_10:end),vq2(idx_max_10:end))) + abs(trapz(v_total_2(1:left_2),vq2(1:left_2))) - abs(trapz(v_total_2(right_2:idx_max_10),vq2(right_2:idx_max_10)));
Wnet_10 = top_area_10 + bottom_area_10;    % [J] Wnet_out for 8 deg
exp_10 = (Wnet_10/Qin_ideal_10)*100;             % [%] Actual efficiency for 8 deg

top_area_12 = abs(trapz(v_total_3(idx_min_12:right_3),vq3(idx_min_12:right_3))) - abs(trapz(v_total_3(left_3:idx_min_12),vq3(left_3:idx_min_12)));
bottom_area_12 = abs(trapz(v_total_3(idx_max_12:end),vq3(idx_max_12:end))) + abs(trapz(v_total_3(1:left_3),vq3(1:left_3))) - abs(trapz(v_total_3(right_3:idx_max_12),vq3(right_3:idx_max_12)));
Wnet_12 = top_area_12 + bottom_area_12;    % [J] Wnet_out for 8 deg
exp_12 = (Wnet_12/Qin_ideal_12)*100;             % [%] Actual efficiency for 8 deg

Wnet_out = [Wnet_8, Wnet_10, Wnet_12]';
exp = [exp_8, exp_10, exp_12]';

%% Experimentally Ideal Efficiency (%)
ideal_8 = 100*(1/(1-T_L_8/T_H_8) + 5/(2*log10(Vmax_8/Vmin_8)))^-1;
ideal_10 = 100*(1/(1-T_L_10/T_H_10) + 5/(2*log10(Vmax_10/Vmin_10)))^-1;
ideal_12 = 100*(1/(1-T_L_12/T_H_12) + 5/(2*log10(Vmax_12/Vmin_12)))^-1;

ideal = [ideal_8, ideal_10, ideal_12]';

% %% Win, Wout, Qin Calculation
% Mass of Air 
m1 = Pmax1*Vmin_8/(R*T_L_8);
m2 = Pmax2*Vmin_10/(R*T_L_10); 
m3 = Pmax3*Vmin_12/(R*T_L_12); 

%% Win, Wout, Qin Calculation
Win_8 = m1*R*T_L_8*log10(Vmin_8/Vmax_8);
Win_10 = m2*R*T_L_10*log10(Vmin_10/Vmax_10);
Win_12 = m3*R*T_L_12*log10(Vmin_12/Vmax_12);

Wout_8 = m1*R*T_H_8*log10(Vmax_8/Vmin_8);
Wout_10 = m2*R*T_H_10*log10(Vmax_10/Vmin_10);
Wout_12 = m3*R*T_H_12*log10(Vmax_12/Vmin_12);

Win = -[Win_8, Win_10, Win_12]';
Wout = [Wout_8, Wout_10, Wout_12]';
Qin = [Qin_ideal_8, Qin_ideal_10, Qin_ideal_12]';
Qout = Qin - Wnet_out; 


%% Energy & Efficiency Table 
fprintf('    Table 1: Efficiency Comparison\n')
colName1 = ["Carnot (%)" "Ideal (%)" "Experimental (%)"];
rowName1 = ["Data 1 ( 8 deg)" "Data 2 (10 deg)" "Data 3 (12 deg)"];
table1 = array2table([carnot, ideal, exp],"RowNames", rowName1, "VariableNames", colName1);
disp(table1)

fprintf('    Table 2: Energy Comparison (x10^-4 except Qin & Qout)\n')
colName1 = ["Win (J)" "Wout (J)" "Wnet_out (J)" "Qin (J)" "Qout (J)"];
rowName1 = ["Data 1 ( 8 deg)" "Data 2 (10 deg)" "Data 3 (12 deg)"];
table1 = array2table([Win*10^4, Wout*10^4, Wnet_out*10^4, Qin, Qout],"RowNames", rowName1, "VariableNames", colName1);
disp(table1)

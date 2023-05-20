%% ASEN 3113 Design Lab
% Blake Wilson, Dakota Harris, Liam Mazzota, Hide Nakanishi, Nathan Evans
% Last Modified: 05/04/2023

%Husekeeping
clc
clear all;
close all;

% Define constants
sigma = 5.67*10^(-8);
alpha_initial = 0.2;        % Absorptivity
epsilon_initial = 0.85;     % Emissivity

% IR backload [W/m^2]
q_ir_backload_winter = 88;
q_ir_backload_summer = 63;
q_ir_backload_eclipse = 11;
q_ir_average = 0.5*(q_ir_backload_summer + q_ir_backload_winter);

% Distances from sun [m]
r_summer = 152.1e9;
r_winter = 147.1e9;
r_equinox = (r_summer + r_winter) *.5;

% Luminosity of sun adjusted for distance [W/m^2]
L_summer = ((3.9e26))/(4*pi*r_summer^2);
L_winter = ((3.9e26))/(4*pi*r_winter^2);
L_equinox = (3.9e26)/(4*pi*r_equinox^2);

% Max affective Area
a_effective_summer = cosd(-23.5);
a_effective_winter = cosd(23.5);
ts = 30 + 273.15; %k
t_surr = 2.7;%k

%% Area Calculations 
% Heat generated by the instrument
interm = (epsilon_initial * sigma *(ts^4 - (t_surr)^4));

%% For summer 
Area_summer_on = (20)/((interm) - (alpha_initial * L_summer * a_effective_summer) - (epsilon_initial * q_ir_backload_summer));


%% For winter
Area_winter_on = (20)/((interm) - (alpha_initial * L_winter * a_effective_winter) - (epsilon_initial * q_ir_backload_winter));


%% For equinox
Area_equinox_on = (20)/((interm) - (alpha_initial * L_equinox * 1) - (epsilon_initial * q_ir_average));

% A_max = 0.289 m^2

%% Modeling solar radiation flux 
ts = 25 + 273.15; %k
interm = (epsilon_initial * sigma *(ts^4 - (t_surr)^4));

Area = Area_winter_on;

% Vector of theta values
theta_plot = linspace(0,(2*pi),1000);
%Vector of time values
time = linspace(0,(24*3600),1000);

eclipse_angle1 = deg2rad(8.685);
eclipse_angle2 = deg2rad(360-8.685);

% Loop to determine solar radiation flux 
for i = 1:length(theta_plot)
    if eclipse_angle1 <= theta_plot(i) && theta_plot(i) <= pi 
        Q_solar_rad_flux_e(i) = (L_equinox) .* alpha_initial .* cosd(0).* sin(theta_plot(i));
    elseif pi < theta_plot(i) && theta_plot(i) <= eclipse_angle2
        Q_solar_rad_flux_e(i) = -(L_equinox) .* alpha_initial .* cosd(0).* sin(theta_plot(i));
    else 
        Q_solar_rad_flux_e(i) = 0;
    end
    
    if theta_plot(i) <= pi
        Q_solar_rad_flux(i) = (L_summer) .* alpha_initial .* cosd(-23.5) .* sin(theta_plot(i)); % FOR SUMMER !!! 
        Q_solar_rad_flux_w(i) = (L_winter) .* alpha_initial .* cosd(23.5) .* sin(theta_plot(i));
    else
        Q_solar_rad_flux(i) = -(L_summer) .* alpha_initial .* cosd(-23.5) .* sin(theta_plot(i)); % FOR SUMMER !!! 
        Q_solar_rad_flux_w(i) = -(L_winter) .* alpha_initial .* cosd(23.5) .* sin(theta_plot(i));
    end

end 

figure()
hold on
grid minor
plot(time/3600,Q_solar_rad_flux,'Color', [0.6350 0.0780 0.1840])
plot(time/3600,Q_solar_rad_flux_w, 'Color',[0 0.4470 0.7410])
plot(time/3600,Q_solar_rad_flux_e, 'Color', [0.4660 0.6740 0.1880])
xlabel('Time (Hours)')
ylabel('Solar Radiation Flux [W/m^2]')
title('Solar Radiation Flux vs. Time during a Day')
xticks([0 6 12 18 24])
xticklabels({'0', '6', '12', '18', '24'})
xlim([0 24])
xline([6 12 18 24])
yline(0)
legend('Summer','Winter','Equinox')
hold off

%% Unheated Temperature & Heater Power Summer

% Heated Power Summer
for idx = 1:length(time)
    if 6*3600 <= time(idx) && time(idx) < 18*3600
        Ts_summer_unheated(idx) = ((20/Area + Q_solar_rad_flux(idx) +q_ir_backload_summer .* epsilon_initial)./(epsilon_initial*sigma) + t_surr^4).^(1/4)-273.15;
        if Ts_summer_unheated(idx) >= 20 && Ts_summer_unheated(idx)<= 30
            Summer_Heater_Power(idx) = 0;
        else
            % q_dot reqired to be in the range
            ts = 20 +273.15;
            interm = (epsilon_initial * sigma *(ts^4 - (t_surr)^4));
            Summer_Heater_Power(idx) = (Area * (interm - Q_solar_rad_flux(idx) - (q_ir_backload_summer * epsilon_initial)))-20;
        end
        
    else
        Ts_summer_unheated(idx) = ((Q_solar_rad_flux(idx) +q_ir_backload_summer .* epsilon_initial)./(epsilon_initial*sigma) + t_surr^4).^(1/4)-273.15; 
        if Ts_summer_unheated(idx) >= -40
            Summer_Heater_Power(idx) = 0;
        else
            % q_dot reqired to be in the range
            ts = -40 + 273.15;
            interm = (epsilon_initial * sigma *(ts^4 - (t_surr)^4));
            Summer_Heater_Power(idx) = Area * (interm - Q_solar_rad_flux(idx) - (q_ir_backload_summer * epsilon_initial));
        end
    end
end

Qtotal_Heater_summer = trapz(time,Summer_Heater_Power);

figure()
hold on
grid minor
title("Summer Temperature and Power Required during a Day")
yyaxis left
p11 = plot(time/3600,Ts_summer_unheated,'r');
ylabel("Unheated Summer Temperature [\circ C]")
yline([30 20 -40],'--',{'Max (Operational)','Min (Operational)','Min (Survival)'},'Color','g')
ylim([-100 40])
yline(0)
ax = gca;
ax.YColor = 'r';
xlabel("Time [Hrs]")


yyaxis right
ylabel("Heater Power [W]")
p12 = plot(time/3600, Summer_Heater_Power,'b');
ylim([-20 70])
yline(0)
ax.YColor = 'b';
xticks([0 6 12 18 24])
xticklabels({'0', '6', '12', '18', '24'})
xlim([0 24])
xline([6 12 18 24])

legend([p11,p12],["Temperature", "Power"])
legend('Location', 'Southeast')
hold off


%% Unheated Temperature & Heater Power Winter Solstice

% Heated Power Winter
for idx = 1:length(time)
    if 6*3600 <= time(idx) && time(idx) < 18*3600
        Ts_winter_unheated(idx) = ((20/Area + Q_solar_rad_flux(idx) +q_ir_backload_winter .* epsilon_initial)./(epsilon_initial*sigma) + t_surr^4).^(1/4)-273.15;
        if Ts_winter_unheated(idx) >= 20 && Ts_winter_unheated(idx)<= 30
            Winter_Heater_Power(idx) = 0;
        else
            % q_dot reqired to be in the range
            ts = 20 +273.15;
            interm = (epsilon_initial * sigma *(ts^4 - (t_surr)^4));
            Winter_Heater_Power(idx) = (Area * (interm - Q_solar_rad_flux(idx) - (q_ir_backload_winter * epsilon_initial)))-20;
        end
        
    else
        Ts_winter_unheated(idx) = ((Q_solar_rad_flux(idx) +q_ir_backload_winter .* epsilon_initial)./(epsilon_initial*sigma) + t_surr^4).^(1/4)-273.15; 
        if Ts_winter_unheated(idx) >= -40
            Winter_Heater_Power(idx) = 0;
        else
            % q_dot reqired to be in the range
            ts = -40 + 273.15;
            interm = (epsilon_initial * sigma *(ts^4 - (t_surr)^4));
            Winter_Heater_Power(idx) = Area * (interm - Q_solar_rad_flux(idx) - (q_ir_backload_winter * epsilon_initial));
        end
    end
end

Qtotal_Heater_winter = trapz(time,Winter_Heater_Power);

figure()
hold on
grid minor
title("Winter Temperature and Power Required during a Day")
yyaxis left
p21 = plot(time/3600,Ts_winter_unheated,'r');
ylabel("Unheated Winter Temperature [\circ C]")
yline([30 20 -40],'--',{'Max (Operational)','Min (Operational)','Min (Survival)'},'Color','g')
ylim([-100 40])
yline(0)
ax = gca;
ax.YColor = 'r';
xlabel("Time [Hrs]")

yyaxis right
ylabel("Heater Power [W]")
p22 = plot(time/3600, Winter_Heater_Power,'b');
ylim([-20 70])
yline(0)
ax.YColor = 'b';
xticks([0 6 12 18 24])
xticklabels({'0', '6', '12', '18', '24'})
xlim([0 24])
xline([6 12 18 24])

legend([p21,p22],["Temperature", "Power"])
legend('Location', 'Southeast')
hold off

%% Unheated Temperature & Heater Power Equinox
% Heated Power Equinox
for idx = 1:length(time)
    if 6*3600 <= time(idx) && time(idx) < 18*3600
        Ts_equinox_unheated(idx) = ((20/Area + Q_solar_rad_flux(idx) +q_ir_average .* epsilon_initial)./(epsilon_initial*sigma) + t_surr^4).^(1/4)-273.15;
        if Ts_equinox_unheated(idx) >= 20 && Ts_equinox_unheated(idx)<= 30
            Equinox_Heater_Power(idx) = 0;
        else
            % q_dot reqired to be in the range
            ts = 20 +273.15;
            interm = (epsilon_initial * sigma *(ts^4 - (t_surr)^4));
            Equinox_Heater_Power(idx) = (Area * (interm - Q_solar_rad_flux(idx) - (q_ir_average * epsilon_initial)))-20;
        end
        
    else
        Ts_equinox_unheated(idx) = ((Q_solar_rad_flux(idx) +q_ir_average .* epsilon_initial)./(epsilon_initial*sigma) + t_surr^4).^(1/4)-273.15; 
        if Ts_equinox_unheated(idx) >= -40
            Equinox_Heater_Power(idx) = 0;
        else
            % q_dot reqired to be in the range
            ts = -40 + 273.15;
            interm = (epsilon_initial * sigma *(ts^4 - (t_surr)^4));
            Equinox_Heater_Power(idx) = Area * (interm - Q_solar_rad_flux(idx) - (q_ir_average * epsilon_initial));
        end
    end
end

Qtotal_Heater_equinox = trapz(time,Equinox_Heater_Power);

figure()
hold on
grid minor
title("Equinox Temperature and Power Required during a Day")
yyaxis left
p31 = plot(time/3600,Ts_equinox_unheated,'r');
ylabel("Unheated Equinox Temperature [\circ C]")
yline([30 20 -40],'--',{'Max (Operational)','Min (Operational)','Min (Survival)'},'Color','g')
ylim([-100 40])
yline(0)
ax = gca;
ax.YColor = 'r';
xlabel("Time [Hrs]")

yyaxis right
ylabel("Heater Power [W]")
p32 = plot(time/3600, Equinox_Heater_Power,'b');
ylim([-20 70])
yline(0)
ax.YColor = 'b';
xticks([0 6 12 18 24])
xticklabels({'0', '6', '12', '18', '24'})
xlim([0 24])
xline([6 12 18 24])

legend([p31,p32],["Temperature", "Power"])
legend('Location', 'Southeast')
hold off

%% New Coating 
% (Aclar Film (Aluminum Backing) 2 mil)
% epsilon = .62;
% alpha = .11;

%Mylar Film ( 3.0 mil aluminum backing) %0.76; 0.183
epsilon = 0.82;
alpha = 0.188;

% New Area Required 
ts_new = 30 + 273.15;
interm_new = (epsilon * sigma *(ts_new^4 - (t_surr)^4));
New_Area = (20)/((interm_new) - (alpha * L_winter * a_effective_winter) - (epsilon * q_ir_backload_winter));

% Loop to determine solar radiation flux 
for i = 1:length(theta_plot)
    if eclipse_angle1 <= theta_plot(i) && theta_plot(i) <= pi 
        Q_solar_rad_flux_e_new(i) = (L_equinox) .* alpha .* cosd(0).* sin(theta_plot(i));
    elseif pi < theta_plot(i) && theta_plot(i) <= eclipse_angle2
        Q_solar_rad_flux_e_new(i) = -(L_equinox) .* alpha .* cosd(0).* sin(theta_plot(i));
    else 
        Q_solar_rad_flux_e_new(i) = 0;
    end
    
    if theta_plot(i) <= pi
        Q_solar_rad_flux_new(i) = (L_summer) .* alpha .* cosd(-23.5) .* sin(theta_plot(i)); % FOR SUMMER !!! 
        Q_solar_rad_flux_w_new(i) = (L_winter) .* alpha .* cosd(23.5) .* sin(theta_plot(i));
    else
        Q_solar_rad_flux_new(i) = -(L_summer) .* alpha .* cosd(-23.5) .* sin(theta_plot(i)); % FOR SUMMER !!! 
        Q_solar_rad_flux_w_new(i) = -(L_winter) .* alpha .* cosd(23.5) .* sin(theta_plot(i));
    end

end 

figure()
hold on
grid minor
plot(time/3600,Q_solar_rad_flux_new,'Color', [0.6350 0.0780 0.1840],'LineWidth',1)
plot(time/3600,Q_solar_rad_flux, '--','Color', [0.6350 0.0780 0.1840])
plot(time/3600,Q_solar_rad_flux_w_new, 'Color',[0 0.4470 0.7410],'LineWidth',1)
plot(time/3600,Q_solar_rad_flux_w, '--', 'Color',[0 0.4470 0.7410])
plot(time/3600,Q_solar_rad_flux_e_new, 'Color', [0.4660 0.6740 0.1880],'LineWidth',1)
plot(time/3600,Q_solar_rad_flux_e, '--', 'Color', [0.4660 0.6740 0.1880])
xlabel('Time [Hrs]')
ylabel('Solar Radiation Flux [W/m^2]')
title('Solar Radiation Flux vs. Time during a Day (Improved Coating)')
xticks([0 6 12 18 24])
xticklabels({'0', '6', '12', '18', '24'})
xlim([0 24])
xline([6 12 18 24])
yline(0)
legend('Summer (Improved)','Summer (Initial)','Winter (Improved)','Winter (Initial)', 'Equinox (Improved)', 'Equinox (Initial)')
legend('Location', 'north')
hold off


%% NEW Unheated Temperature for the operational mode at Summer - Heater power for the operational mode
% Heated Power Summer
for idx = 1:length(time)
    if 6*3600 <= time(idx) && time(idx) < 18*3600
        Ts_summer_unheated_new(idx) = ((20/New_Area + Q_solar_rad_flux(idx) +q_ir_backload_summer .* epsilon)./(epsilon*sigma) + t_surr^4).^(1/4)-273.15;
        if Ts_summer_unheated_new(idx) >= 20 && Ts_summer_unheated_new(idx)<= 30
            Summer_Heater_Power_new(idx) = 0;
        else
            % q_dot reqired to be in the range
            ts = 20 +273.15;
            interm = (epsilon * sigma *(ts^4 - (t_surr)^4));
            Summer_Heater_Power_new(idx) = (New_Area * (interm - Q_solar_rad_flux(idx) - (q_ir_backload_summer * epsilon)))-20;
        end
        
    else
        Ts_summer_unheated_new(idx) = ((Q_solar_rad_flux(idx) +q_ir_backload_summer .* epsilon)./(epsilon*sigma) + t_surr^4).^(1/4)-273.15; 
        if Ts_summer_unheated_new(idx) >= -40
            Summer_Heater_Power_new(idx) = 0;
        else
            % q_dot reqired to be in the range
            ts = -40 + 273.15;
            interm = (epsilon * sigma *(ts^4 - (t_surr)^4));
            Summer_Heater_Power_new(idx) = New_Area * (interm - Q_solar_rad_flux(idx) - (q_ir_backload_summer * epsilon));
        end
    end
end

Qtotal_Heater_summer_new = trapz(time,Summer_Heater_Power_new);
Qtotal_summer = [Qtotal_Heater_summer; Qtotal_Heater_summer_new].*10^-6;

figure()
hold on
grid minor
title("Summer Temperature and Power Required during a Day (Improved Coating)")
yyaxis left
p41 = plot(time/3600,Ts_summer_unheated_new,'r','LineWidth',1);
p42 = plot(time/3600,Ts_summer_unheated,'--r');
ylabel("Unheated Summer Temperature [\circ C]")
yline([30 20 -40],'--',{'Max (Operational)','Min (Operational)','Min (Survival)'},'Color','g')
ylim([-100 40])
yline(0)
ax = gca;
ax.YColor = 'r';
xlabel("Time [Hrs]")

yyaxis right
ylabel("Heater Power [W]")
p43 = plot(time/3600, Summer_Heater_Power_new,'b','LineWidth',1);
p44 = plot(time/3600, Summer_Heater_Power,'--b');
ylim([-20 70])
yline(0)
ax.YColor = 'b';
xticks([0 6 12 18 24])
xticklabels({'0', '6', '12', '18', '24'})
xlim([0 24])
xline([6 12 18 24])

legend([p41,p42,p43,p44],["Temperature (Improved)" "Temperature (Initial)", "Power (Improved)", "Power (Initial)"])
legend('Location', 'Southeast')
hold off


%% NEW Unheated Temperature for the operational mode at Winter - Heater power for the operational mode
% Heated Power Winter
for idx = 1:length(time)
    if 6*3600 <= time(idx) && time(idx) < 18*3600
        Ts_winter_unheated_new(idx) = ((20/New_Area + Q_solar_rad_flux(idx) +q_ir_backload_winter .* epsilon)./(epsilon*sigma) + t_surr^4).^(1/4)-273.15;
        if Ts_winter_unheated_new(idx) >= 20 && Ts_winter_unheated_new(idx)<= 30
            Winter_Heater_Power_new(idx) = 0;
        else
            % q_dot reqired to be in the range
            ts = 20 +273.15;
            interm = (epsilon * sigma *(ts^4 - (t_surr)^4));
            Winter_Heater_Power_new(idx) = (New_Area * (interm - Q_solar_rad_flux(idx) - (q_ir_backload_winter * epsilon)))-20;
        end
        
    else
        Ts_winter_unheated_new(idx) = ((Q_solar_rad_flux(idx) +q_ir_backload_winter .* epsilon)./(epsilon*sigma) + t_surr^4).^(1/4)-273.15; 
        if Ts_winter_unheated_new(idx) >= -40
            Winter_Heater_Power_new(idx) = 0;
        else
            % q_dot reqired to be in the range
            ts = -40 + 273.15;
            interm = (epsilon * sigma *(ts^4 - (t_surr)^4));
            Winter_Heater_Power_new(idx) = New_Area * (interm - Q_solar_rad_flux(idx) - (q_ir_backload_winter * epsilon));
        end
    end
end

Qtotal_Heater_winter_new = trapz(time,Winter_Heater_Power_new);
Qtotal_winter = [Qtotal_Heater_winter; Qtotal_Heater_winter_new].*10^-6;

figure()
hold on
grid minor
title("Winter Temperature and Power Required during a Day (Improved Coating)")
yyaxis left
p51 = plot(time/3600,Ts_winter_unheated_new,'r','LineWidth',1);
p52 = plot(time/3600,Ts_winter_unheated,'--r');
ylabel("Unheated Winter Temperature [\circ C]")
yline([30 20 -40],'--',{'Max (Operational)','Min (Operational)','Min (Survival)'},'Color','g')
ylim([-100 40])
yline(0)
ax = gca;
ax.YColor = 'r';
xlabel("Time [Hrs]")

yyaxis right
ylabel("Heater Power [W]")
p53 = plot(time/3600, Winter_Heater_Power_new,'b','LineWidth',1);
p54 = plot(time/3600, Winter_Heater_Power,'--b');
ylim([-20 70])
yline(0)
ax.YColor = 'b';
xticks([0 6 12 18 24])
xticklabels({'0', '6', '12', '18', '24'})
xlim([0 24])
xline([6 12 18 24])

legend([p51,p52,p53,p54],["Temperature (Improved)" "Temperature (Initial)", "Power (Improved)", "Power (Initial)"])
legend('Location', 'Southeast')
hold off

%% NEW Unheated Temperature for the operational mode at Equinox - Heater power for the operational mode
% Heated Power Equinox
for idx = 1:length(time)
    if 6*3600 <= time(idx) && time(idx) < 18*3600
        Ts_equinox_unheated_new(idx) = ((20/New_Area + Q_solar_rad_flux(idx) +q_ir_average .* epsilon)./(epsilon*sigma) + t_surr^4).^(1/4)-273.15;
        if Ts_equinox_unheated_new(idx) >= 20 && Ts_equinox_unheated_new(idx)<= 30
            Equinox_Heater_Power_new(idx) = 0;
        else
            % q_dot reqired to be in the range
            ts = 20 +273.15;
            interm = (epsilon * sigma *(ts^4 - (t_surr)^4));
            Equinox_Heater_Power_new(idx) = (New_Area * (interm - Q_solar_rad_flux(idx) - (q_ir_average * epsilon)))-20;
        end
        
    else
        Ts_equinox_unheated_new(idx) = ((Q_solar_rad_flux(idx) +q_ir_average .* epsilon)./(epsilon*sigma) + t_surr^4).^(1/4)-273.15; 
        if Ts_equinox_unheated_new(idx) >= -40
            Equinox_Heater_Power_new(idx) = 0;
        else
            % q_dot reqired to be in the range
            ts = -40 + 273.15;
            interm = (epsilon * sigma *(ts^4 - (t_surr)^4));
            Equinox_Heater_Power_new(idx) = New_Area * (interm - Q_solar_rad_flux(idx) - (q_ir_average * epsilon));
        end
    end
end

Qtotal_Heater_equinox_new = trapz(time,Equinox_Heater_Power_new);
Qtotal_equinox = [Qtotal_Heater_equinox; Qtotal_Heater_equinox_new].*10^-6;

figure()
hold on
grid minor
title("Equinox Temperature and Power Required during a Day (Improved Coating)")
yyaxis left
p61 = plot(time/3600,Ts_equinox_unheated_new,'r','LineWidth',1);
p62 = plot(time/3600,Ts_equinox_unheated,'--r');
ylabel("Unheated Equinox Temperature [\circ C]")
yline([30 20 -40],'--',{'Max (Operational)','Min (Operational)','Min (Survival)'},'Color','g')
ylim([-100 40])
yline(0)
ax = gca;
ax.YColor = 'r';
xlabel("Time [Hrs]")

yyaxis right
ylabel("Heater Power [W]")
p63 = plot(time/3600, Equinox_Heater_Power_new,'b','LineWidth',1);
p64 = plot(time/3600, Equinox_Heater_Power,'--b');
ylim([-20 70])
yline(0)
ax.YColor = 'b';
xticks([0 6 12 18 24])
xticklabels({'0', '6', '12', '18', '24'})
xlim([0 24])
xline([6 12 18 24])

legend([p61,p62,p63,p64],["Temperature (Improved)" "Temperature (Initial)", "Power (Improved)", "Power (Initial)"])
legend('Location', 'Southeast')
hold off


%% Temp with heater on
% % Summer solctice Temp with heating
% Ts_summer_heated = (((20 + Summer_Heater_Power_Operational)/Area + Q_solar_rad_flux +q_ir_backload_summer .* epsilon_initial)./(epsilon_initial*sigma) + t_surr^4).^(1/4)-273.15;
% Ts_summer_heated_survival = ((Summer_Heater_Power_Survival/Area + Q_solar_rad_flux +q_ir_backload_summer .* epsilon_initial)./(epsilon_initial*sigma) + t_surr^4).^(1/4)-273.15; 
% 
% % Winter solctice Temp with heating
% Ts_winter_heated_Operational = (((20+Winter_Heater_Power_Operational)/Area + Q_solar_rad_flux_w +q_ir_backload_winter .* epsilon_initial)./(epsilon_initial*sigma) + t_surr^4).^(1/4)-273.15;
% Ts_winter_heated_survival = ((Winter_Heater_Power_Survival/Area + Q_solar_rad_flux_w +q_ir_backload_winter .* epsilon_initial)./(epsilon_initial*sigma) + t_surr^4).^(1/4)-273.15; 
% 
% % Equinox Temp with heating
% Ts_equinox_heated_Operational = (((20+Equinox_Heater_Power_Operational)/Area + Q_solar_rad_flux_e + q_ir_average .* epsilon_initial)./(epsilon_initial*sigma) + t_surr^4).^(1/4)-273.15;
% Ts_equinox_heated_survival = ((Equinox_Heater_Power_Survival/Area + Q_solar_rad_flux_e +q_ir_average.* epsilon_initial)./(epsilon_initial*sigma) + t_surr^4).^(1/4)-273.15; 
% 
% % Equinox Temp with heating and new coating
% Ts_equinox_heated_Operational_new = (((20+Equinox_Heater_Power_Operational_new)/New_Area_winter_on + Q_solar_rad_flux_e_new + q_ir_average .* epsilon)./(epsilon*sigma) + t_surr^4).^(1/4)-273.15;
% Ts_equinox_heated_survival_new = ((Equinox_Heater_Power_Survival_new/New_Area_winter_on + Q_solar_rad_flux_e_new +q_ir_average.* epsilon)./(epsilon*sigma) + t_surr^4).^(1/4)-273.15;
% 
%     for i = 1:length(Ts_equinox_heated_survival_new)
%     if Ts_equinox_heated_survival_new(i) < -40
%         Ts_equinox_heated_survival_new(i) = -40;
%     end
%     end



% figure()
% hold on
% grid minor
% plot(time/3600, Ts_summer_heated_Operational)
% plot(time/3600,Ts_winter_heated_Operational)
% plot(time/3600,Ts_equinox_heated_Operational)
% plot(time/3600,Ts_equinox_heated_Operational_new)
% xlabel("Time [Hrs]")
% ylabel("Temperature [\circ C]")
% title("Operational: Heated Temperature vs Time during a Day")
% yline([30 20],'--',{'Max (Operational)','Min (Operational)'},'Color','g')
% xticks([0 6 12 18 24])
% xticklabels({'6 am', '12 pm', '6 pm', '0 am', '6 am'})
% xlim([0 24])
% ylim([19 31])
% xline([6 12 18 24])
% legend("Summer Solstace","Winter Solstace","Equinox","Equinox: New Coating")
% 
% 
% figure ()
% hold on
% grid minor
% plot(time/3600,Ts_summer_heated_survival)
% plot(time/3600,Ts_winter_heated_survival)
% plot(time/3600,Ts_equinox_heated_survival)
% plot(time/3600,Ts_equinox_heated_survival_new)
% xlabel("Time [Hrs]")
% ylabel("Temperature [\circ C]")
% legend("Summer Solstace","Winter Solstace","Equinox","Equinox: New Coating")
% yline(-40,'--',{'Min (Survival)'},'Color','g')
% 
% title("Survival: Heated Temperature vs Time during a Day")
% xticks([0 6 12 18 24])
% xticklabels({'6 am', '12 pm', '6 pm', '0 am', '6 am'})
% xlim([0 24])
% xline([6 12 18 24])
% legend("Summer Solstace","Winter Solstace","Equinox","Equinox: New Coating")


%% Table 
% Table 1: Area Calculation
fprintf('Table 1: Minimum Satellite Area Required\n')
colName1 = [" Initial " " Improved "];
rowName1 = ["Area [m^2]"];
table1 = array2table([Area, New_Area],"RowNames", rowName1, "VariableNames", colName1);
disp(table1)

% Table 2: Total Power Required by the Heater for the Operational Mode 
Qtotal = [Qtotal_summer(1)+Qtotal_winter(1)+Qtotal_equinox(1); Qtotal_summer(2)+Qtotal_winter(2)+Qtotal_equinox(2)];
fprintf('Table 2: Total Energy Required by Heater\n')
colName1 = [" Summer " " Winter " " Equinox " " Total "];
rowName1 = ["Original [MJ]", "Improved [MJ]"];
table1 = array2table([Qtotal_summer, Qtotal_winter, Qtotal_equinox, Qtotal],"RowNames", rowName1, "VariableNames", colName1);
disp(table1)




%% INITIAL COATING Source
%https://spark.iop.org/suns-luminosity 
% USING WINTER !!! .2891
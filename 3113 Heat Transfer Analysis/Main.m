%% ASEN 3113 Lab 02

clc;
close all;
clear all;

steel_plot = readmatrix('Steel_21V_192mA');
brass26_plot = readmatrix('Brass_26V_245mA');
brass29_plot = readmatrix('Brass_29V_273mA');
aluminum28_plot = readmatrix('Aluminum_28V_270mA');
aluminum30_plot = readmatrix('Aluminum_30V_290mA');


index_steel = process('Steel_21V_192mA', 30,1);
index_brass26 = process('Brass_26V_245mA', 25,1);
index_brass29 = process('Brass_29V_273mA',25,1);
index_aluminum28 = process('Aluminum_28V_270mA',15,0);
index_aluminum30 = process('Aluminum_30V_290mA',20,0);

ss_steel = [steel_plot(index_steel,3); steel_plot(index_steel,4); steel_plot(index_steel,5); steel_plot(index_steel,6); steel_plot(index_steel,7); steel_plot(index_steel,8); steel_plot(index_steel,9); steel_plot(index_steel,10)];
[T0_s, Hexp_s, Han_s] = steadystate(ss_steel,21,192,16.2) % mA converted in func

ss_brass26 = [brass26_plot(index_brass26,3); brass26_plot(index_brass26,4); brass26_plot(index_brass26,5); brass26_plot(index_brass26,6); brass26_plot(index_brass26,7); brass26_plot(index_brass26,8); brass26_plot(index_brass26,9); brass26_plot(index_brass26,10)];
[T0_b26, Hexp_b26, Han_b26] = steadystate(ss_brass26,26,245,115)% mA converted in func

ss_brass29 = [brass29_plot(index_brass29,3); brass29_plot(index_brass29,4); brass29_plot(index_brass29,5); brass29_plot(index_brass29,6); brass29_plot(index_brass29,7); brass29_plot(index_brass29,8); brass29_plot(index_brass29,9); brass29_plot(index_brass29,10)];
[T0_b29, Hexp_b29, Han_b29] = steadystate(ss_brass29,29,273,115) % mA converted in func

ss_aluminum28 = [aluminum28_plot(index_aluminum28,2);aluminum28_plot(index_aluminum28,3);aluminum28_plot(index_aluminum28,4);aluminum28_plot(index_aluminum28,5);aluminum28_plot(index_aluminum28,6);aluminum28_plot(index_aluminum28,7);aluminum28_plot(index_aluminum28,8);aluminum28_plot(index_aluminum28,9);];
[T0_a28, Hexp_a28, Han_a28] = steadystate(ss_aluminum28,28,270,130) % mA converted in func

ss_aluminum30 = [aluminum30_plot(index_aluminum30,2);aluminum30_plot(index_aluminum30,3);aluminum30_plot(index_aluminum30,4);aluminum30_plot(index_aluminum30,5);aluminum30_plot(index_aluminum30,6);aluminum30_plot(index_aluminum30,7);aluminum30_plot(index_aluminum30,8);aluminum30_plot(index_aluminum30,9);];
[T0_a30, Hexp_a30, Han_a30] = steadystate(ss_aluminum30,30,290,130)
% mA converted in func





%% Time Dependant Temperature Profiles
% Recall from the lab discussion that the general solution for ��(��, ��) involves taking the summation of �� from 1 to ∞, and
% that the solution converges and can be approximated very accurately by the summation of �� over a finite range. Generate
% plots of the analytical temperature ��(��, ��) — using the analytical steady state slope ������ — as a function of time for each
% thermocouple location in a single color, overlayed onto the experimental temperature data for each thermocouple
% location in a different color. We will call this first version of the analytical model “Model IA.”

% Determine needed variables
L = 5.875 * .0254; %[m]
xvec = [(1+3/8):.5:4.875] * .0254;
t_s = steel_plot(:,1);
t_b26 = brass26_plot(:,1);
t_b29 = brass29_plot(:,1);
t_a28 = aluminum28_plot(:,1);
t_a30 = aluminum30_plot(:,1);
alpha_s = 4.05e-06;
alpha_b26 = 3.593e-05;
alpha_b29 = 3.593e-05;
alpha_a = 4.819e-05;

% Steel
% model IA
for i = 1:length(xvec)
    
    [TempVec] = TempProfile(Han_s,10,xvec(i),L,T0_s,alpha_s,t_s);
    figure(4)
    hold on
    grid on
    plot(t_s,TempVec,'r')
    plot(t_s,steel_plot(:,2+i),'b')
end
figure(4)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Steel: Experimental vs Analytical resuluts using H analytical')
legend("Analytical Model using Han","Experimental Data")

%Model IB
for i = 1:length(xvec)
    
    [TempVec] = TempProfile(Hexp_s,10,xvec(i),L,T0_s,alpha_s,t_s);
    figure(5)
    hold on
    grid on
    plot(t_s,TempVec,'r')
    plot(t_s,steel_plot(:,2+i),'b')
end
figure(5)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Steel: Experimental vs Analytical resuluts using H exp')
legend("Analytical Model using Hexp","Experimental Data")





% Brass 26
% model IA
for i = 1:length(xvec)
    
    [TempVec] = TempProfile(Han_b26,10,xvec(i),L,T0_b26,alpha_b26,t_b26);
    figure(6)
    hold on
    grid on
    plot(t_b26,TempVec,'r')
    plot(t_b26,brass26_plot(:,2+i),'b')
end
figure(6)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Brass 26V: Experimental vs Analytical resuluts using H analytical')
legend("Analytical Model using Han","Experimental Data")

%Model IB
for i = 1:length(xvec)
    
    [TempVec] = TempProfile(Hexp_b26,10,xvec(i),L,T0_b26,alpha_b26,t_b26);
    figure(7)
    hold on
    grid on
    plot(t_b26,TempVec,'r')
    plot(t_b26,brass26_plot(:,2+i),'b')
end
figure(7)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Brass 26V: Experimental vs Analytical resuluts using H exp')
legend("Analytical Model using Hexp","Experimental Data")





% Brass 29
% model IA
for i = 1:length(xvec)
    
    [TempVec] = TempProfile(Han_b29,10,xvec(i),L,T0_b29,alpha_b29,t_b29);
    figure(8)
    hold on
    grid on
    plot(t_b29,TempVec,'r')
    plot(t_b29,brass29_plot(:,2+i),'b')
end
figure(8)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Brass 29V: Experimental vs Analytical resuluts using H analytical')
legend("Analytical Model using Han","Experimental Data")

%Model IB
for i = 1:length(xvec)
    
    [TempVec] = TempProfile(Hexp_b29,10,xvec(i),L,T0_b29,alpha_b29,t_b29);
    figure(9)
    hold on
    grid on
    plot(t_b29,TempVec,'r')
    plot(t_b29,brass29_plot(:,2+i),'b')
end
figure(9)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Brass 29V: Experimental vs Analytical resuluts using H exp')
legend("Analytical Model using Hexp","Experimental Data")





% Aluminum 28V
% model IA
for i = 1:length(xvec)
    
    [TempVec] = TempProfile(Han_a28,10,xvec(i),L,T0_a28,alpha_a,t_a28);
    figure(10)
    hold on
    grid on
    plot(t_a28,TempVec,'r')
    plot(t_a28,aluminum28_plot(:,1+i),'b')
end
figure(10)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Aluminum 28V: Experimental vs Analytical resuluts using H analytical')
legend("Analytical Model using Han","Experimental Data")

%Model IB
for i = 1:length(xvec)
    
    [TempVec] = TempProfile(Hexp_a28,10,xvec(i),L,T0_a28,alpha_a,t_a28);
    figure(11)
    hold on
    grid on
    plot(t_a28,TempVec,'r')
    plot(t_a28,aluminum28_plot(:,1+i),'b')
end
figure(11)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Aluminum 28V: Experimental vs Analytical resuluts using H exp')
legend("Analytical Model using Hexp","Experimental Data")


% Aluminum 30V
% model IA
for i = 1:length(xvec)
    
    [TempVec] = TempProfile(Han_a30,10,xvec(i),L,T0_a30,alpha_a,t_a30);
    figure(12)
    hold on
    grid on
    plot(t_a30,TempVec,'r')
    plot(t_a30,aluminum30_plot(:,1+i),'b')
end
figure(12)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Aluminum 30V: Experimental vs Analytical resuluts using H analytical')
legend("Analytical Model using Han","Experimental Data")

%Model IB
for i = 1:length(xvec)
    
    [TempVec] = TempProfile(Hexp_a30,10,xvec(i),L,T0_a30,alpha_a,t_a30);
    figure(13)
    hold on
    grid on
    plot(t_a30,TempVec,'r')
    plot(t_a30,aluminum30_plot(:,1+i),'b')
end
figure(13)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Aluminum 30V: Experimental vs Analytical resuluts using H exp')
legend("Analytical Model using Hexp","Experimental Data")






%% Initial State Ditribution Model II

L = 5.875 * .0254; %[m]
xvec = [(1+3/8):.5:4.875] * .0254;
t_s = steel_plot(:,1);
t_b26 = brass26_plot(:,1);
t_b29 = brass29_plot(:,1);
t_a28 = aluminum28_plot(:,1);
t_a30 = aluminum30_plot(:,1);
alpha_s = 4.05e-06;
alpha_b26 = 3.593e-05;
alpha_b29 = 3.593e-05;
alpha_a = 4.819e-05;
M_s = 34.589;
T0_Is = 12.692;
M_b26 = 7.03044;
T0_Ib26 = 11.779;
M_b29 = 5.8118;
T0_Ib29 = 11.764;
T0_Ia28 = 10.214;
M_a28 = -2.0622;
T0_Ia30 = 8.5232;
M_a30 = -2.8121;
% Steel
for i = 1:length(xvec)
    
    [TempVec] = TempProfileII(Hexp_s,10,xvec(i),L,T0_Is,alpha_s,t_s,M_s);
    figure(14)
    hold on
    grid on
    txt_an = ['Thermocouple',num2str(i)];
    plot(t_s,TempVec,'r')
    plot(t_s,steel_plot(:,2+i),'b')
end
figure(14)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Steel: Experimental vs Analytical resuluts using Initial State Distribution')
legend("Analytical Model using Hexp","Experimental Data")

%Brass 26
for i = 1:length(xvec)
    
    [TempVec] = TempProfileII(Hexp_b26,10,xvec(i),L,T0_Ib26,alpha_b26,t_b26,M_b26);
    figure(15)
    hold on
    grid on
    txt_an = ['Thermocouple',num2str(i)];
    plot(t_b26,TempVec,'r')
    plot(t_b26,brass26_plot(:,2+i),'b')
end
figure(15)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Brass 26V: Experimental vs Analytical resuluts using Initial State Distribution')
legend("Analytical Model using Hexp","Experimental Data")

% Brass 29
for i = 1:length(xvec)
    
    [TempVec] = TempProfileII(Hexp_b29,10,xvec(i),L,T0_Ib29,alpha_b29,t_b29,M_b29);
    figure(16)
    hold on
    grid on
    txt_an = ['Thermocouple',num2str(i)];
    plot(t_b29,TempVec,'r')
    plot(t_b29,brass29_plot(:,2+i),'b')
end
figure(16)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Brass 29V: Experimental vs Analytical resuluts using Initial State Distribution')
legend("Analytical Model using Hexp","Experimental Data")


% Aluminum 28
for i = 1:length(xvec)
    
    [TempVec] = TempProfileII(Hexp_a28,10,xvec(i),L,T0_Ia28,alpha_a,t_a28,M_a28);
    figure(17)
    hold on
    grid on
    txt_an = ['Thermocouple',num2str(i)];
    plot(t_a28,TempVec,'r')
    plot(t_a28,aluminum28_plot(:,1+i),'b')
end
figure(17)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Aluminum 28V: Experimental vs Analytical resuluts using Initial State Distribution')
legend("Analytical Model using Hexp","Experimental Data")

%Aluminum30
for i = 1:length(xvec)
    
    [TempVec] = TempProfileII(Hexp_a30,10,xvec(i),L,T0_Ia30,alpha_a,t_a30,M_a30);
    figure(18)
    hold on
    grid on
    txt_an = ['Thermocouple',num2str(i)];
    plot(t_a30,TempVec,'r')
    plot(t_a30,aluminum30_plot(:,1+i),'b')
end
figure(18)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Aluminum 30V: Experimental vs Analytical resuluts using Initial State Distribution')
legend("Analytical Model using Hexp","Experimental Data")



%% Varying alpha

L = 5.875 * .0254; %[m]
xvec = [(1+3/8):.5:4.875] * .0254;
t_s = steel_plot(:,1);
t_b26 = brass26_plot(:,1);
t_b29 = brass29_plot(:,1);
t_a28 = aluminum28_plot(:,1);
t_a30 = aluminum30_plot(:,1);
alpha_s = linspace(3e-6, 5e-06,5);
alpha_b = linspace(3e-5,4e-05,5);
alpha_a = linspace(4e-05,5.5e-5,5);
color = ["r" "g" "y" "m" "k"];
%Steel
for j = 1:length(xvec)
    %for i = 1:length(alpha_s)
    
    [TempVec] = TempProfile(Hexp_s,10,xvec(j),L,T0_s,alpha_s(4),t_s);
    figure(19)
    hold on
    grid on
    txt_an = ['Thermocouple',num2str(i)];
    plot(t_s,TempVec,'color',color(1))
    plot(t_s,steel_plot(:,2+j),'b','linewidth',1.5)
    %end
end
figure(19)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Model III: Steel')
legend("","Experimental Data")
% Magenta is closest Alpha analytical is 4.5 e-6

%Brass
for j = 1:length(xvec)
    %for i = 1:length(alpha_b)
    [TempVec] = TempProfile(Hexp_b29,10,xvec(j),L,T0_b29,alpha_b(1),t_b29);
    figure(20)
    hold on
    grid on
    txt_an = ['Thermocouple',num2str(i)];
    plot(t_b29,TempVec,'color',color(1))
    plot(t_b29,brass29_plot(:,2+j),'b','linewidth',1.5)
    %end
end
figure(20)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Model III: Brass')
legend("","Experimental Data")
% red is closest Alpha analytical is 3 e-5


% Aluminum
for j = 1:length(xvec)
    %for i = 1:length(alpha_a)
    [TempVec] = TempProfile(Hexp_a28,10,xvec(j),L,T0_a28,alpha_a(1),t_a28);
    figure(21)
    hold on
    grid on
    txt_an = ['Thermocouple',num2str(i)];
    plot(t_a28,TempVec,'color',color(1))
    plot(t_a28,aluminum28_plot(:,1+j),'b','linewidth',1.5)
    %end
end
figure(21)
xlabel('Time')
ylabel('Temperature [deg C]')
title('Model III: Aluminum')
legend("","Experimental Data")
% red is closest Alpha analytical is 4 e-5


%% Time to steady state:
alpha_adj_steel = 4.5*10^(-6);
alpha_adj_brass = 3.0 * 10^(-5);
alpha_adj_alum = 4.0 * 10^(-5);

[tss_steel, Fo_steel] = fourier_number(alpha_adj_steel);
[tss_brass, Fo_brass] = fourier_number(alpha_adj_brass);
[tss_alum, Fo_alum] = fourier_number(alpha_adj_alum);

%% Th8 Measured Temperature over Time
Hexp_all = [Hexp_s, Hexp_b26, Hexp_b29, Hexp_a28, Hexp_a30];
M_all = [M_s,M_b26,M_b29,M_b29,M_a28,M_a30];
T0_all = [T0_s, T0_b26, T0_b29, T0_a28, T0_a30];
T0_I_all = [T0_Is, T0_Ib26, T0_Ib29, T0_Ia28, T0_Ia30];
alpha_all = [4.5e-6, 3e-5, 3e-5, 4e-5, 4e-5];   % [m^2/s] Adjusted thermal diffusivity
t_all{1} = t_s;
t_all{2} = t_b26;
t_all{3} = t_b29;
t_all{4} = t_a28;
t_all{5} = t_a30;
true_plot{1} = steel_plot(:,10);
true_plot{2} = brass26_plot(:,10);
true_plot{3} = brass29_plot(:,10);
true_plot{4} = aluminum28_plot(:,9);
true_plot{5} = aluminum30_plot(:,9);
name = ["Steel", "Brass 26V", "Brass 29V", "Aluminum 28V", "Aluminum 30V"];
sys_error = 2;  % [˚C] +/- 2 ˚C systematic error

for i = 1:length(Hexp_all)
    figure(i+21)
    Th8_vecII{i,:} = TempProfileII(Hexp_all(i),10,xvec(8),L,T0_I_all(i),alpha_all(i),t_all{i},M_all(i));
    Th8_vecIII{i,:} = TempProfile(Hexp_all(i),10,xvec(8),L,T0_all(i),alpha_all(i),t_all{i});
    plot(t_all{i}, true_plot{i},'LineWidth', 2); hold on
    plot(t_all{i}, Th8_vecII{i},'LineWidth', 2); hold on
    plot(t_all{i}, Th8_vecIII{i},'LineWidth', 2); hold on
    plot(t_all{i}, true_plot{i}+sys_error,'--g','LineWidth', 2); hold on
    plot(t_all{i}, true_plot{i}-sys_error,'--g','LineWidth', 2); hold on
    xlabel('Time (s)')
    ylabel('Temperature at Thermocouple 8 (^\circ C)')
    title(sprintf('%s: Experimental vs Analytical Temperature Change at Thermocouple 8 over Time', name(i)));
    legend('Experimental Model','Analytical Model II', 'Adjusted Analytical Model II', '+/- 2^{\circ}C Systematic Error')
    grid on
end

%% Functions
% This function calculates forier number and tss
function [tss, Fo] = fourier_number(alpha)
fo = 0;
t = 0;
L = 5.875/39.37;

while fo <= .2
    fo = (alpha * t)/(L^2);
    
    
    
    t = t + .1;
end
tss = t;
Fo = fo;
end

% with blips in the data, there are lots of points in the graph in which
% the slope of the line is acctually 0 therefore it is essential to have a
% tolerance and check the other slopes in order to verify the blips aren't
% just being selected due to their 0 slope.
function [index] = process(filename, sub,mat)

data = readmatrix(filename);
pos = [11/8; 15/8; 19/8; 23/8; 27/8; 31/8; 35/8; 39/8];

if mat == 1
    time = data(:,1);
    Tc1 = data(:,3);
    Tc2 = data(:,4);
    Tc3 = data(:,5);
    Tc4 = data(:,6);
    Tc5 = data(:,7);
    Tc6 = data(:,8);
    Tc7 = data(:,9);
    Tc8 = data(:,10);
else
    time = data(:,1);
    Tc1 = data(:,2);
    Tc2 = data(:,3);
    Tc3 = data(:,4);
    Tc4 = data(:,5);
    Tc5 = data(:,6);
    Tc6 = data(:,7);
    Tc7 = data(:,8);
    Tc8 = data(:,9);
end



tolerance = .0001;
index = 1;
for i = (sub+1):length(Tc1)
    x = (Tc1(i) - Tc1(i-sub))/(time(i) - time(i-sub));
    x2 = (Tc2(i) - Tc2(i-sub))/(time(i) - time(i-sub));
    x3 = (Tc3(i) - Tc3(i-sub))/(time(i) - time(i-sub));
    x4 = (Tc4(i) - Tc4(i-sub))/(time(i) - time(i-sub));
    x5 = (Tc5(i) - Tc5(i-sub))/(time(i) - time(i-sub));
    x6 = (Tc6(i) - Tc6(i-sub))/(time(i) - time(i-sub));
    x7 = (Tc7(i) - Tc7(i-sub))/(time(i) - time(i-sub));
    x8 = (Tc8(i) - Tc8(i-sub))/(time(i) - time(i-sub));
    
    
    if (x == 0)...
            &&  (-tolerance <= x2) && x2 <= tolerance &&  (-tolerance <= x3) && x3 <= tolerance && (-tolerance <= x4) && x4 <= tolerance && (-tolerance <= x5)...
            && x5 <= tolerance && (-tolerance <= x6) && x6 <= tolerance && (-tolerance <= x7) && x7 <= tolerance && (-tolerance <= x8) && x8 <= tolerance
        index = i;
        break;
    end
end
end


% in order to compute the steady state slope, we will pass in 1 array with
% all of our temperatures at the steady state index. Then using matalab
% line interp we will find the slope and y intercept. Rn we only care about
% the slope and this is the H we will be finding.

%array must be passed in Th0 -> Th8
function [T0, hexp, han] = steadystate(array,V,A,K)
pos = [11/8; 15/8; 19/8; 23/8; 27/8; 31/8; 35/8; 39/8]./39.37;

inter = polyfit(pos,array,1);
hexp = inter(1); % degrees C/m
T0 = inter(2); %degrees c

area = pi* (0.0127^2); %m^2
q_dot = (V * A)/1000; % W

%k is W/ (m * K)

han = q_dot/(area * K);
end



function [TempVec] = TempProfile(H,n,x,L,T0,alpha,t)
% Creates time dependant temp profile where outputs are a Temperature vector
% over time t and vector of time t.
% Inouts:
%       H: Steady state slope
%       n: number of iterations
%       x: location of thermo couple
%       L: length of bar
%       T0: temp at cool end of bar
%       t: time range
    b_n = @(n,H) -(8*H*L*(-1)^(n+1))/(2*n*pi-pi)^2;
    lambda_n = @(n) ((2*n-1)*pi)/(2*L);


    sum =0;
    for j = 1:n %loop for n
        sum = sum + b_n(j,H) .* sin(lambda_n(j).*x).*exp(-lambda_n(j).^2.*alpha.*t);

        TempVec = T0 + H.*x + sum;
    end

end



function [TempVecII] = TempProfileII(H,n,x,L,T0,alpha,t,M)
    b_n = @(n,H) (8*(M-H)*L*(-1)^(n+1))/(2*n*pi-pi)^2;
    lambda_n = @(n) ((2*n-1)*pi)/(2*L);


    sum =0;
    for j = 1:n %loop for n
        sum = sum + b_n(j,H) .* sin(lambda_n(j).*x).*exp(-lambda_n(j).^2.*alpha.*t);

        TempVecII = T0 + H.*x + sum;
    end
end


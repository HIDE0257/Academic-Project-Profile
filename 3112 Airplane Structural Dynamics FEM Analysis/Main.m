clear all; clc; close all;
Constants2 = getConstants2();

L = Constants2.L; 
% Assembling master mass matrix
M2_p1 = (Constants2.cM2.*[19272 1458*Constants2.L 5928 -642*Constants2.L 0 0;
                         1458*Constants2.L 172*Constants2.L^2 642*Constants2.L -73*Constants2.L^2 0 0;
                         5928 642*Constants2.L 38544 0 5928 -642*Constants2.L;
                         -642*Constants2.L -73*Constants2.L^2 0 344*Constants2.L^2 642*Constants2.L -73*Constants2.L^2;
                         0 0 5928 642*Constants2.L 19272 -1458*Constants2.L;
                         0 0 -642*Constants2.L -73*Constants2.L^2 -1458*Constants2.L 172*Constants2.L^2]); % First part of Master mass matrix
M2_p2 = [0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 Constants2.Mt Constants2.St;
      0 0 0 0 Constants2.St Constants2.It]; % Second part of Master mass matrix for 2 element model
M2 = M2_p1 + M2_p2; % Full master mass equation for 2 element model
% Assembling Master stiffness equation
K2 = Constants2.cK2.*[24 6*Constants2.L -24 6*Constants2.L 0 0;
                      6*Constants2.L 2*Constants2.L^2 -6*Constants2.L Constants2.L^2 0 0 
                      -24 -6*Constants2.L 48 0 -24 6*Constants2.L
                      6*Constants2.L Constants2.L^2 0 4*Constants2.L^2 -6*Constants2.L Constants2.L^2
                      0 0 -24 -6*Constants2.L 24 -6*Constants2.L
                      0 0 6*Constants2.L Constants2.L^2 -6*Constants2.L 2*Constants2.L^2]; % Master stiffness matrix for 2 element model
% Reduced Mass and stiffness equations:
M2_reduced = M2((3:end),(3:end)); % Reduced master mass equation for two element model
K2_reduced = K2((3:end),(3:end)); % Reduced master stiffness equation for two element model
[eigvec2,eig2] = eig(K2_reduced,M2_reduced); % Eigenvalues/vectors of modeled dynamical system

% Resizing eigenvector to 6x6:
zeromatrix = zeros(6,6);
zeromatrix(3:6,3:6) = eigvec2;
eigvec2 = zeromatrix;

M4_p1 = Constants2.cM4.*[77088 2916*L 23712 -1284*L 0 0 0 0 0 0; 2916*L 172*L^2 11284*L -73*L^2 0 0 0 0 0 0; 
    23712 1284*L 154176 0 23712 -1284*L 0 0 0 0; 
    -1284*L -73*L^2 0 344*L^2 1284*L -73*L^2 0 0 0 0;
    0 0 23712 1284*L 154176 0 23712 -1284*L 0 0; 
    0 0 -1284*L -73*L^2 0 344*L^2 1284*L -73*L^2 0 0; 
    0 0 0 0 23712 1284*L 154176 0 23712 -1284*L; 
    0 0 0 0 -1284*L -73*L^2 0 344*L^2 1284*L -73*L^2; 
    0 0 0 0 0 0 23712 1284*L 77099 -2916*L; 
    0 0 0 0 0 0 -1284*L -73*L^2 -2916*L 172*L^2];
M4_p2 = [0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 Constants2.Mt Constants2.St;
         0 0 0 0 0 0 0 0 Constants2.St Constants2.It]; % Second part of master mass equation
M4 = M4_p1 + M4_p2; % Full Master mass equation for four element model
% Assembling Master stiffness matrix:
K4 = Constants2.cK4.*[96 12*Constants2.L -96 12*Constants2.L 0 0 0 0 0 0;
                      12*Constants2.L 2*(Constants2.L^2) -12*Constants2.L Constants2.L^2 0 0 0 0 0 0;
                      -96 -12*Constants2.L 192 0 -96 12*Constants2.L 0 0 0 0;
                      12*Constants2.L Constants2.L^2 0 4*Constants2.L^2 -12*Constants2.L Constants2.L^2 0 0 0 0;
                      0 0 -96 -12*Constants2.L 192 0 -96 12*Constants2.L 0 0;
                      0 0 12*Constants2.L Constants2.L^2 0 4*Constants2.L^2 -12*Constants2.L Constants2.L^2 0 0;
                      0 0 0 0 -96 -12*Constants2.L 192 0 -96 12*Constants2.L;
                      0 0 0 0 12*Constants2.L Constants2.L^2 0 4*Constants2.L^2 -12*Constants2.L Constants2.L^2;
                      0 0 0 0 0 0 -96 -12*Constants2.L 96 -12*Constants2.L;
                      0 0 0 0 0 0 12*Constants2.L Constants2.L^2 -12*Constants2.L 2*Constants2.L^2]; % Master stiffness matrix for four element model
% Creating Eigenvalue problem:
M4_reduced = M4((3:end),(3:end)); % Reduced master mass equation for four element model
K4_reduced = K4((3:end),(3:end)); % Reduced stiffness equation for four element model
[eigvec4,eig4] = eig(K4_reduced,M4_reduced); % Eigenvalues/vectors of modeled dynamical system
% Re-sizing eigenvector with 0's (8x8)
zeromatrix = zeros(10,10);
zeromatrix(3:10,3:10) = eigvec4;
eigvec4 = zeromatrix;

%%
filename = ["test_2min_all_8" "test_2min_nose_2" "test_2min_tail_1" "test_2min_wing_1" "test_5min_all_2"];
loaction = ["All_2" "All_5" "Wing" "Nose" "Tail"];
%% All Points of intrest
% 2 minute test
peakAll_5 = plotExperimentalData(filename(1),loaction(1));
% 5 minute test
peakAll_2 =plotExperimentalData(filename(5),loaction(2));
eigenAll = sort(mean([peakAll_2;peakAll_5]));
eigenValTail_V = eigenAll;

peakTail_H = plotExperimentalData(filename(3),loaction(5));
peakTail_H = mean(peakTail_H);
eigval = [eigenValTail_V(1:2),peakTail_H];
w_tailV = 2*pi*eigval;

A1_tailV = K2_reduced - w_tailV(1)^2.*M2_reduced;
A2_tailV = K2_reduced - w_tailV(2)^2.*M2_reduced;
A3_tailV = K2_reduced - w_tailV(3)^2.*M2_reduced;
[V1_tailV,~] = eig(A1_tailV);
[V2_tailV,~] = eig(A2_tailV);
[V3_tailV,~] = eig(A3_tailV);
zeromatrix11 = zeros(6,6);
zeromatrix11(3:6,3:6) = V1_tailV;
V1_tailV = zeromatrix11;
zeromatrix22 = zeros(6,6);
zeromatrix22(3:6,3:6) = V2_tailV;
V2_tailV = zeromatrix22;
zeromatrix33 = zeros(6,6);
zeromatrix33(3:6,3:6) = V3_tailV;
V3_tailV = zeromatrix33;

a1_tailV = K4_reduced - w_tailV(1)^2.*M4_reduced;
a2_tailV = K4_reduced - w_tailV(2)^2.*M4_reduced;
a3_tailV = K4_reduced - w_tailV(3)^2.*M4_reduced;
[v1_tailV,~] = eig(a1_tailV);
[v2_tailV,~] = eig(a2_tailV);
[v3_tailV,~] = eig(a3_tailV);
zeromatrix1 = zeros(10,10);
zeromatrix1(3:10,3:10) = v1_tailV;
v1_tailV = zeromatrix1;
zeromatrix2 = zeros(10,10);
zeromatrix2(3:10,3:10) = v2_tailV;
v2_tailV = zeromatrix2;
zeromatrix3 = zeros(10,10);
zeromatrix3(3:10,3:10) = v3_tailV;
v3_tailV = zeromatrix3;

% ~ 12Hz resonance for vertiacl tail
% ~ 24Hz resonance for wing
% ~ 42Hz resonance for Nose
% ~ 46Hz resonance for Horisontal tail


% Defining inputs for ploteigenvector function:
%% Two-Element Model
L2 = Constants2.L; % length of FEM Mode
ne2 = 2; % Number of beam finite elements in the model
nsub2 = 10; % Number of subdivisions
scale2 = 1; % Scaling value for eigenvectors

% Experiemntals
[X1_e,V1_e] = ploteigenvector(L2,V1_tailV(:,3),ne2,nsub2,scale2,16); % Plot for eigenvector 1
[X2_e,V2_e] = ploteigenvector(L2,V2_tailV(:,4),ne2,nsub2,scale2,17); % Plot for eigenvector 1
[X3_e,V3_e] = ploteigenvector(L2,V3_tailV(:,5),ne2,nsub2,scale2,18); % Plot for eigenvector 1
V1_e = V1_e/norm(V1_e(:,1));
V2_e = V2_e/norm(V2_e(:,1));
V3_e = V3_e/norm(V3_e(:,1));

% FEM
[X1,V1] = ploteigenvector(L2,eigvec2(:,3),ne2,nsub2,scale2,1); % Plot for eigenvector 1
[X2,V2] = ploteigenvector(L2,eigvec2(:,4),ne2,nsub2,scale2,2); % Plot for eigenvector 2
[X3,V3] = ploteigenvector(L2,eigvec2(:,5),ne2,nsub2,scale2,3); % Plot for eigenvector 3
V1 = V1/norm(V1(:,1));
V2 = V2/norm(V2(:,1));
V3 = V3/norm(V3(:,1));

dx1_E2 = V1(2) - V1(1);
dx2_E2 = V2(2) - V2(1);
dx3_E2 = V3(2) - V3(1);

%% Four-Element Model 
L4 = Constants2.L; % length of FEM Mode
ne4 = 4; % Number of beam finite elements in the model
nsub4 = 10; % Number of subdivisions
scale4 = 1; % Scaling value for eigenvectors
% Plotting eigenvector for four element model:
% Experiemntals
[x1_e,v1_e] = ploteigenvector(L4,v1_tailV(:,3),ne4,nsub4,scale4,16); % Plot for eigenvector 1
[x2_e,v2_e] = ploteigenvector(L4,v2_tailV(:,4),ne4,nsub4,scale4,17); % Plot for eigenvector 1
[x3_e,v3_e] = ploteigenvector(L4,v3_tailV(:,5),ne4,nsub4,scale4,18); % Plot for eigenvector 1
v1_e = v1_e/norm(v1_e(:,1));
v2_e = v2_e/norm(v2_e(:,1));
v3_e = v3_e/norm(v3_e(:,1));

% FEM
[x1,v1] = ploteigenvector(L4,eigvec4(:,3),ne4,nsub4,scale4,1); % Plot for eigenvector 1
[x2,v2] = ploteigenvector(L4,eigvec4(:,4),ne4,nsub4,scale4,2); % Plot for eigenvector 2
[x3,v3] = ploteigenvector(L4,eigvec4(:,5),ne4,nsub4,scale4,3); % Plot for eigenvector 3
v1 = v1/norm(v1(:,1));
v2 = v2/norm(v2(:,1));
v3 = v3/norm(v3(:,1));

dx1_E4 = v1(2) - v1(1);
dx2_E4 = v2(2) - v2(1);
dx3_E4 = v3(2) - v3(1);

% [~,E2_v2_e] = ploteigenvector(L2,V3_tailV(:,4),ne2,nsub2,scale2,18); 
% E2_v2_e = E2_v2_e/norm(E2_v2_e(:,1),2);
% 
% [~,E2_v2] = ploteigenvector(L2,eigvec2(:,4),ne2,nsub2,scale2,18); 
% E2_v2 = E2_v2/norm(E2_v2(:,1),2);
% 
% [~,E4_v2_e] = ploteigenvector(L4,v3_tailV(:,4),ne4,nsub4,scale4,18);
% E4_v2_e = E4_v2_e/norm(E4_v2_e(:,1),2);
% [~,E4_v2] = ploteigenvector(L4,eigvec4(:,4),ne4,nsub4,scale4,18); 
% E4_v2 = E4_v2/norm(E4_v2(:,1),2);

%% Error Calculations
% Error of FEM and Experiment for the 2 Element
error1 = 100*(abs(trapz(X1(:,1),V1(:,1))) - abs(trapz(X1(:,1),V1_e(:,1))))/abs(trapz(X1(:,1),V1(:,1)));
error2 = 100*(abs(trapz(X2(:,1),V2(:,1))) - abs(trapz(X2(:,1),V2_e(:,1))))/abs(trapz(X2(:,1),V2(:,1)));
error3 = 100*(abs(trapz(X3(:,1),V3(:,1))) - abs(trapz(X3(:,1),V3_e(:,1))))/abs(trapz(X3(:,1),V3(:,1)));

% error11 = sum(abs((V1(:,1) - V1_e(:,1))))*dx1_E2*100;
% error22 = sum(abs((V2(:,1) - V2_e(:,1))))*dx2_E2*100;
% error33 = sum(abs((V3(:,1) - V3_e(:,1))))*dx3_E2*100;

% Error of FEM and Experiment for the 4 Element
error11 = 100*(abs(trapz(x1(:,1),v1(:,1))) - abs(trapz(x1(:,1),v1_e(:,1))))/abs(trapz(x1(:,1),v1(:,1)));
error22 = 100*(abs(trapz(x2(:,1),v2(:,1))) - abs(trapz(x2(:,1),v2_e(:,1))))/abs(trapz(x2(:,1),v2(:,1)));
error33 = 100*(abs(trapz(x3(:,1),v3(:,1))) - abs(trapz(x3(:,1),v3_e(:,1))))/abs(trapz(x3(:,1),v3(:,1)));

% error1 = sum(abs((v1(:,1) - v1_e(:,1))))*dx1_E4*100;
% error2 = sum(abs((v2(:,1) - v2_e(:,1))))*dx2_E4*100;
% error3 = sum(abs((v3(:,1) - v3_e(:,1))))*dx3_E4*100;

e_E2 = [error1; error2; error3];
e_E4 = [error11; error22; error33];

colName2 = [ "2 Elements (%)" "4 Elements (%)"];
rowName2 = ["Mode Shape 1" "Mode Shape 2" "Mode Shape 3"];
Q2_table = array2table([e_E2,e_E4],"RowNames", rowName2, "VariableNames", colName2);
disp(Q2_table)

%% Plots 
figure(1)
hold on
plt1 = plot(x1(:,1),-v1(:,1),'LineWidth',1.5,'Color', 'r');
plt1_e = plot(x1_e(:,1),v1_e(:,1),'LineWidth',1.5,'Color', 'b');
title('Four-Element Model: Mode Shapes')
xlabel('Length of Mode [inches]')
ylabel('Shape');
legend([plt1(1),plt1_e(1)],'FEM Mode Shape 1', 'Experimental Mode Shape 1')
hold off

figure(2)
hold on
plt2 = plot(x2(:,1),-v2(:,1),'LineWidth',1.5,'Color','r');
plt2_e = plot(x2_e(:,1),v2_e(:,1),'LineWidth',1.5,'Color', 'b');
title('Four-Element Model: Mode Shapes')
xlabel('Length of Mode [inches]')
ylabel('Shape')
legend([plt2(1),plt2_e(1)],'FEM Mode Shape 2', 'Experimental Mode Shape 2')
hold off

figure(3)
hold on
plt2 = plot(x3(:,1),-v3(:,1),'LineWidth',1.5,'Color','r');
plt2_e = plot(x3_e(:,1),-v3_e(:,1),'LineWidth',1.5,'Color', 'b');
title('Four-Element Model: Mode Shapes')
xlabel('Length of Mode [inches]')
ylabel('Shape')
legend([plt2(1),plt2_e(1)],'FEM Mode Shape 3', 'Experimental Mode Shape 3')
hold off

figure(4)
hold on 
plot(x1_e(:,1),-v1_e(:,1),'LineWidth',1.5,'Color','r');
title('Four-Element Model: Experimental Mode Shape 1')
xlabel('Length of Mode [inches]')
ylabel('Shape')
hold off

figure(5)
hold on 
plot(x2_e(:,1),-v2_e(:,1),'LineWidth',1.5,'Color','r');
title('Four-Element Model: Experimental Mode Shape 2')
xlabel('Length of Mode [inches]')
ylabel('Shape')
hold off

figure(6)
hold on 
plot(x3_e(:,1),-v3_e(:,1),'LineWidth',1.5,'Color','r');
title('Four-Element Model: Experimental Mode Shape 3')
xlabel('Length of Mode [inches]')
ylabel('Shape')
hold off

% 2 Elements
figure(10)
hold on
plt11 = plot(X1(:,1),V1(:,1),'LineWidth',1.5,'Color', 'r');
plt11_e = plot(X1_e(:,1),V1_e(:,1),'LineWidth',1.5,'Color', 'b');
title('Two-Element Model: Mode Shapes')
xlabel('Length of Mode [inches]')
ylabel('Shape');
legend([plt11(1),plt11_e(1)],'FEM Mode Shape 1', 'Experimental Mode Shape 1')
hold off

figure(11)
hold on
plt22 = plot(X2(:,1),V2(:,1),'LineWidth',1.5,'Color','r');
plt22_e = plot(X2_e(:,1),V2_e(:,1),'LineWidth',1.5,'Color', 'b');
title('Two-Element Model: Mode Shapes')
xlabel('Length of Mode [inches]')
ylabel('Shape')
legend([plt22(1),plt22_e(1)],'FEM Mode Shape 2', 'Experimental Mode Shape 2')
hold off

figure(12)
hold on
plt11 = plot(X3(:,1),V3(:,1),'LineWidth',1.5,'Color', 'r');
plt11_e = plot(X3_e(:,1),V3_e(:,1),'LineWidth',1.5,'Color', 'b');
title('Two-Element Model: Mode Shapes')
xlabel('Length of Mode [inches]')
ylabel('Shape');
legend([plt11(1),plt11_e(1)],'FEM Mode Shape 3', 'Experimental Mode Shape 3')
hold off

figure(13)
hold on 
plot(X1_e(:,1),V1_e(:,1),'LineWidth',1.5,'Color','r');
title('Two-Element Model: Experimental Mode Shape 1')
xlabel('Length of Mode [inches]')
ylabel('Shape')
hold off

figure(14)
hold on 
plot(X2_e(:,1),V2_e(:,1),'LineWidth',1.5,'Color','r');
title('Two-Element Model: Experimental Mode Shape 2')
xlabel('Length of Mode [inches]')
ylabel('Shape')
hold off

figure(15)
hold on 
plot(X3_e(:,1),V3_e(:,1),'LineWidth',1.5,'Color','r');
title('Two-Element Model: Experimental Mode Shape 3')
xlabel('Length of Mode [inches]')
ylabel('Shape')
hold off

figure(20)
P1 = plot(X1_e(:,1),V1_e(:,1),'LineWidth',1.5,'Color','r'); hold on
p1 = plot(x1_e(:,1),v1_e(:,1),'LineWidth',1.5,'Color','b'); hold on
title('Theoretical Mode Shape 3')
xlabel('Length of Mode [inches]')
ylabel('Shape')
legend([P1(1),p1(1)],'Two-Element Model', 'Four-Element Model')

figure(21)
P2 = plot(X2_e(:,1),V2_e(:,1),'LineWidth',1.5,'Color','r'); hold on
p2 = plot(x2_e(:,1),v2_e(:,1),'LineWidth',1.5,'Color','b'); hold on
title('Theoretical Mode Shape 3')
xlabel('Length of Mode [inches]')
ylabel('Shape')
legend([P2(1),p2(1)],'Two-Element Model', 'Four-Element Model')

figure(22)
P3 = plot(X3_e(:,1),V3_e(:,1),'LineWidth',1.5,'Color','r'); hold on
p3 = plot(x3_e(:,1),-v3_e(:,1),'LineWidth',1.5,'Color','b'); hold on
title('Theoretical Mode Shape 3')
xlabel('Length of Mode [inches]')
ylabel('Shape')
legend([P3(1),p3(1)],'Two-Element Model', 'Four-Element Model')


classdef slp
    properties (Constant)
        e = 0.7;
        e0 = 0.7;
        C_fe = 0.003;
    end
    methods (Static)
        function [slope_lift,CL_R] = calc_slope(S_wet,S_ref,b)
            %flat plate data
            FP_Langle = [-13.0728 ,-12.09 ,-11.0735 ,-10.0574,-9.07559 ,-8.09424 ,-7.14725 ,-6.13293,-5.05147 ,-4.07163 ,-3.09146 ,-2.07747,0.930398 ,1.97819 ,2.99231 ,4.00636 ,5.05408 ,7.0479 ,8.02842 ,9.04358 ,10.0595 ,11.0079 ,12.0582,13.0068];
            FP_Cl = [-0.787245,-0.804155,-0.815366,-0.806654,-0.783716,-0.740854,-0.678074,-0.598204,-0.492707,-0.387227,-0.295978,-0.201877,0.0889578,0.185911,0.274319,0.365573,0.465372,0.664953,0.74197,0.787684,0.802089,0.805097,0.799584,0.794054];
            
            FP_Dangle = [-14.035,-13.0991,-12.0962,-11.0932,-10.09,-8.08331,-7.11346,-5.17358,-1.12963,2.04325,4.94718,5.91449,6.94837,7.9823,8.98279,9.95014,10.9509,11.9853,12.9862,13.9871,15.9898,17.9919,23.2319];
            FP_Cd = [0.201482,0.192235, 0.179017, 0.164474, 0.145957, 0.106276, 0.0877583, 0.0480747 ,0.0104336, 0.0125044,0.0403958, 0.0596267, 0.0821705, 0.104052, 0.126595, 0.145164, 0.163071, 0.17833, 0.194251, 0.210833, 0.230092, 0.259284 ,0.315714];
            
            FP_Cd_aoa = interp1(FP_Dangle, FP_Cd, FP_Langle);
            
            %BEst fit line for Cl vs aoa flat plate
            Cl_SlopeLine = polyfit(FP_Langle(5:19),FP_Cl(5:19),1); %polyfit output slope and y intercept
            Cl_Slope = @(x) Cl_SlopeLine(1)*x + Cl_SlopeLine(2);%y= mx+b
            
            
            %calculating 3d wing Cl
            AR = b^2/S_ref;
            slope_lift = Cl_SlopeLine(1)/(1+((57.3*Cl_SlopeLine(1))/(pi*slp.e*AR))); %finding a
            

            alpha0 = (-1)*Cl_SlopeLine(2)/Cl_SlopeLine(1);%finding alpha at 0
            alpha = FP_Langle;
            wing_3dCl = slope_lift * (alpha - alpha0);%final 3d wing Cl calculation
            
            %3d wing cd
            CD_wing = FP_Cd_aoa + ((wing_3dCl.^2)/(pi*slp.e*AR));
            
            %aircraft cd
            CD0 = slp.C_fe*(S_wet/S_ref);
            %e_o = 1.78*(1-(0.045*(AR^0.68)))-0.64;%raymers method
            k = 1/(pi*slp.e0*AR);
            CL_R = sqrt(CD0/k);           % Maximum CL : CD0 = kCL^2
            
            cd3dmin = min(CD_wing);
            p = find(CD_wing==cd3dmin);
            Cl_air = wing_3dCl(p);
            
            CD_aircraft = CD0+ k*((wing_3dCl - Cl_air).^2);
            
%             %PLOT aoa vs lift
%             figure(3)
%             plot(FP_Langle, wing_3dCl)
%             hold on
%             legend('2d truth', 'best fit', '3d flat plate', 'Horizontal Tail Lift Curve');
%             xlabel('angle of attack');
%             ylabel('Coefficient of Lift');
%             grid on
            
%             %PLOT DRAG POLAR 3D
%             figure(4)
%             %3d wing plot
%             plot (FP_Cl,FP_Cd_aoa)
%             hold on
%             plot (wing_3dCl, CD_wing)
%             hold on
%             plot (wing_3dCl, CD_aircraft)
%             hold on
%             
%             ylabel('Coefficient of Drag');
%             xlabel('Coefficient of Lift');
%             legend('2d Flat Plate Wing','3d Flat Plate Wing','3D Aircraft', '3D Flat Plate Tail', '3D Tail')
%             grid on
         end
    end
end

% figure(3)
% plot(FP_Langle, clcd)
% plot(FP_Dangle, FP_Cd)
% xlabel('angle of attack');
% ylabel('Coefficient of Drag');
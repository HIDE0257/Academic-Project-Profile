function [c_l, c_dw] = DiamondAirfoil(M,alpha,epsilon1, epsilon2)
gamma = 1.4;    % Thermally + Calorically perfect gas
N = length(alpha);
c_l = NaN(N,1);
c_dw = NaN(N,1);
for i = 1:N
    if 0 <= abs(alpha(i)) && abs(alpha(i)) < epsilon1   % 0 <= |alpha| < epsilon1
        % Define geometric angles
        if 0 <= alpha(i)    % When angle of attack is POSITIVE with respect to freestream
            theta2 = epsilon1 - alpha(i);   % [deg] Anlge of the left upper surface
            theta3 = epsilon1 + epsilon2;   % [deg] Angle of the right upper surface
            theta4 = epsilon1 + alpha(i);   % [deg] Angle of the left lower surface
            theta5 = theta3;                % [deg] Angle of the right lower surface
        else    % When angle of attack is NEGATIVE
            theta2 = epsilon1 + abs(alpha(i));  % [deg] Anlge of the left upper surface
            theta3 = epsilon1 + epsilon2;       % [deg] Angle of the right upper surface
            theta4 = epsilon1 - abs(alpha(i));  % [deg] Angle of the left lower surface
            theta5 = theta3;                    % [deg] Angle of the right lower surface
        end
        
        % Region 1-2 & 1-4 (Oblique Shocks)
        beta2 = ObliqueShockBeta(M,theta2,gamma,'Weak');    % [deg] Shock angle (upper) at the leading edge
        beta4 = ObliqueShockBeta(M,theta4,gamma,'Weak');    % [deg] Shock angle (lower) at the leading edge
        M12_n = M*sind(beta2);  % Mach normal before shock at region 1-2
        M14_n = M*sind(beta4);  % Mach normal before shock at region 1-4
        % Get pressure ratio (P(local)/P1) and Mach normal after shock
        % [P_24] = [P2/P1, P4/P1]   [M24_n] = [M2_n, M4_n]
        [~,~,P_24,~,M24_n,~,~] = flownormalshock(gamma,[M12_n,M14_n],'mach');
        % [M_24] = [M2, M4]
        M_24 = M24_n./[sind(beta2 - theta2),sind(beta4 - theta4)];  % Convert normal to M along the streamline at the regions
        
        % Region 2-3 & 4-5 (Prandtl-Meyer Expanssion Fans)
        % [nu1_35] = [nu1_3, nu1_5]
        [~,nu1_35,~] = flowprandtlmeyer(gamma,M_24,'mach'); % [deg] Prandtl-Meyer angle at region 3 & 5
        nu2_35 = [theta3,theta5] + nu1_35;  % [deg] New Prandtl_Meyer angle after the fans
        [M_35(1),~,~] = flowprandtlmeyer(gamma,nu2_35(1),'nu'); % M3
        [M_35(2),~,~] = flowprandtlmeyer(gamma,nu2_35(2),'nu'); % M5
        % [P_2435] = [P2/P02, P4/P04, P3/P03, P5/P05]
        [~,~,P_2435,~,~] = flowisentropic(gamma,[M_24,M_35],'mach');
        
        % Final pressure ratios
        P2P1 = P_24(1);                     % P2/P1
        P3P1 = P_2435(3)*P2P1/P_2435(1);    % P3/P1
        P4P1 = P_24(2);                     % P4/P1
        P5P1 = P_2435(4)*P4P1/P_2435(2);    % P4/P1
        
    elseif epsilon1 == abs(alpha(i))    % |alpha| = epsilon1
        if 0 < alpha(i)    % POSITIVE angle
            theta2 = 0;                     % [deg] Anlge of the left upper surface
            theta3 = epsilon1 + epsilon2;   % [deg] Anlge of the right upper surface
            theta4 = alpha(i) + epsilon1;   % [deg] Anlge of the left lower surface
            theta5 = theta3;                % [deg] Anlge of the right lower surface
            
            % Region 1-4 (Oblique Shock)
            beta = ObliqueShockBeta(M,theta4,gamma,'Weak');
            theta = theta4;
        else    % NEGATIVE angle
            theta2 = abs(alpha(i)) + epsilon1;  % [deg] Anlge of the left upper surface
            theta3 = epsilon1 + epsilon2;       % [deg] Anlge of the right upper surface
            theta4 = 0;                         % [deg] Anlge of the left lower surface
            theta5 = theta3;                    % [deg] Anlge of the right lower surface
            
            % Region 1-2 (Oblique Shock)
            beta = ObliqueShockBeta(M,theta2,gamma,'Weak');
            theta = theta2;
        end
        
        M1_n = M*sind(beta);    % M1 normal
        % Pressure ratio and Mach: [P2/P1, M2 normal] or [P4/P1, M4 normal]
        [~,~,P_2or4,~,M2or4_n,~,~] = flownormalshock(gamma,M1_n,'mach');
        M2or4 = M2or4_n/sind(beta - theta); % Mach: [M2] or [M4]
        
        % Region 1-2 or 1-4 & 2-3 & 4-5 (Prandtl-Meyer Expanssion Fans)
        if theta == theta4  % If POSITIVE angle
            % µ1 = PM(gamma,[M1, M4])
            [~,nu1_35,~] = flowprandtlmeyer(gamma,[M,M2or4],'mach');
        elseif theta == theta2 % If NEGATIVE angle
            % µ1 = PM(gamma, [M2, M1])
            [~,nu1_35,~] = flowprandtlmeyer(gamma,[M2or4,M],'mach');
        end
        
        nu2_35 = [theta3,theta5] + nu1_35;  % µ2 = [theta3, theta5] + [µ1_3, µ1_5]
        [M_35(1),~,~] = flowprandtlmeyer(gamma,nu2_35(1),'nu'); % M3
        [M_35(2),~,~] = flowprandtlmeyer(gamma,nu2_35(2),'nu'); % M5
        
        if theta == theta4  % If POSITIVE angle
            % Pressure ratio = [P1/P01, P4/P04, P3/P03, P5/P05] =
            % f(gamma,[M1,M4,M3,M5])
            [~,~,P_1435,~,~] = flowisentropic(gamma,[M,M2or4,M_35],'mach');
            P2P1 = 0;                       % P2/P1
            P3P1 = P_1435(3)/P_1435(1);     % P3/P1
            P4P1 = P_2or4;                  % P4/P1
            P5P1 = P_1435(4)*P4P1/P_1435(2);% P5/P1
        elseif theta == theta2  % If POSITIVE angle
            % Pressure ratio = [P2/P02, P1/P01, P3/P03, P5/P03] =
            % f(gamma,[M2,M1,M3,M5]) *P4/P04 = P1/P01
            [~,~,P_2135,~,~] = flowisentropic(gamma,[M2or4,M,M_35],'mach');
            P2P1 = P_2or4;                  % P2/p1
            P3P1 = P_2135(3)*P2P1/P_2135(1);% P3/P1
            P4P1 = 0;                       % P4/P1
            P5P1 = P_2135(4)/P_2135(2);     % P5/P1
        end
        
    else    % |alpha| > epsilon1
        if 0 < alpha(i)    % If POSITIVE angle
            theta2 = alpha(i) - epsilon1;   % [deg] Anlge of the left upper surface
            theta3 = epsilon1 + epsilon2;   % [deg] Anlge of the right upper surface
            theta4 = alpha(i) + epsilon1;   % [deg] Anlge of the left lower surface
            theta5 = theta3;                % [deg] Anlge of the right lower surface

            % Region 1-4 (Oblique Shock)
            beta = ObliqueShockBeta(M,theta4,gamma,'Weak');
            if beta >= 0 && isreal(beta)  % Ignore non-real, negative numbers
                M1_n = M*sind(beta);
                [~,~,P_4,~,M4_n,~,~] = flownormalshock(gamma,M1_n,'mach');
                M4 = M4_n/sind(beta - theta4);
                if M4 >= 1  % Mach4 has to be >= 1
                    % Region 1-2 & 2-3 & 4-5 (Prandtl-Meyer Expanssion Fans)
                    [~,nu1_1,~] = flowprandtlmeyer(gamma,M,'mach');
                    nu1_2 = theta2 + nu1_1;
                    [M2,~,~] = flowprandtlmeyer(gamma,nu1_2,'nu');
                    [~,nu1_35,~] = flowprandtlmeyer(gamma,[M2,M4],'mach');
                    nu2_35 = [theta3,theta5] + nu1_35;
                    [M3,~,~] = flowprandtlmeyer(gamma,nu2_35(1),'nu');
                    [M5,~,~] = flowprandtlmeyer(gamma,nu2_35(2),'nu');
                    % Pressure ratios: [P1/P01, P2/P02, P3/P03, P4/P04, P5/P05] =
                    % f(gamma, [M1,M2,M3,M4,M5])
                    [~,~,P_12345,~,~] = flowisentropic(gamma,[M,M2,M3,M4,M5],'mach');
                    
                    % Final Pressure ratios
                    P2P1 = P_12345(2)/P_12345(1);       % P2/P1
                    P3P1 = P_12345(3)/P_12345(1);       % P3/P1
                    P4P1 = P_4;                         % P4/P1
                    P5P1 = P_12345(5)*P4P1/P_12345(4);  % P5/P1
                else    % M4 < 1 --> Returns to NaN and skip the loop
                    continue
                end
                
            else    % Beta is non-real, negative number --> Returns to NaN and skip the loop
                continue
            end
            
        else    % When angle of attack is NEGATIVE
            theta2 = abs(alpha(i)) + epsilon1;  % [deg] Anlge of the left upper surface
            theta3 = epsilon1 + epsilon2;       % [deg] Anlge of the right upper surface
            theta4 = abs(alpha(i)) - epsilon1;  % [deg] Anlge of the left lower surface
            theta5 = theta3;                    % [deg] Anlge of the right lower surface
            
            % Region 1-2 (Oblique Shock)
            beta = ObliqueShockBeta(M,theta2,gamma,'Weak');
            if beta >= 0 && isreal(beta)  % Ignore non-real, negative numbers
                M1_n = M*sind(beta);
                
                [~,~,P_2,~,M2_n,~,~] = flownormalshock(gamma,M1_n,'mach');
                M2 = M2_n/sind(beta - theta2);
                if M2 >= 1 % M2 has be >= 1
                    % Region 2-3 & 1-4 & 4-5 (Prandtl-Meyer Expassion Fans)
                    [~,nu1_1,~] = flowprandtlmeyer(gamma,M,'mach');
                    nu1_2 = theta4 + nu1_1;
                    [M4,~,~] = flowprandtlmeyer(gamma,nu1_2,'nu');
                    [~,nu1_35,~] = flowprandtlmeyer(gamma,[M2,M4],'mach');
                    nu2_35 = [theta3,theta5] + nu1_35;
                    [M3,~,~] = flowprandtlmeyer(gamma,nu2_35(1),'nu');
                    [M5,~,~] = flowprandtlmeyer(gamma,nu2_35(2),'nu');
                    % Pressure ratios: [P1/P01, P2/P02, P3/P03, P4/P04, P5/P05] =
                    % f(gamma, [M1,M2,M3,M4,M5])
                    [~,~,P_12345,~,~] = flowisentropic(gamma,[M,M2,M3,M4,M5],'mach');
                    
                    P2P1 = P_2;                         % P2/P1
                    P3P1 = P_12345(3)*P2P1/P_12345(2);  % P3/P1
                    P4P1 = P_12345(4)/P_12345(1);       % P4/P1
                    P5P1 = P_12345(5)/P_12345(1);       % P5/P1
                else    % M2 < 1 --> Returns to NaN and skip the loop
                    continue
                end
                
            else    % Beta is non-real, negative number --> Returns to Nan and skip the loop
                continue
            end
        end
    end
    
    Cn_d = gamma*M^2*(tand(epsilon1) + tand(epsilon2)); % Denomerator of lift normal coeff
    Cn_n = 2*((P4P1 - P2P1)*tand(epsilon2) + (P5P1 - P3P1)*tand(epsilon1)); % Numerator
    Cn = Cn_n/Cn_d; % Lift normal coefficient
    
    An_v = tand(epsilon1)*tand(epsilon2)/(tand(epsilon1) + tand(epsilon2));
    An_P = P2P1 + P4P1 - P3P1 - P5P1;
    Ca = (2/(gamma*M^2))*An_P*An_v; % Drag axial coefficient
    
    c_l(i) = Cn*cosd(alpha(i)) - Ca*sind(alpha(i));     % Sectional lift coefficient
    c_dw(i) = Cn*sind(alpha(i)) + Ca*cosd(alpha(i));    % Sectional wave drag coeffcient
end

end

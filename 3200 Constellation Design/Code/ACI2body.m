%% ACItoBody
function r_body = ACI2body(r_ACI,T_benu,t)
% This function rotates a vector in the ACI frame to the body frame

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% r_ACI - Position in ACI frame (Astroid Centerd Inertial Frame)
% T_benu - Peroid of rotation for Benu (seconds_)
% OUTPUT 
%  r_body- Position in body frame

Omega = 2*pi/T_benu; % rad/sec

    for i = 1:length(t)

        theta3 = Omega*t(i); % Rotation angle about benu z-axis

        % Rotation Matrix
        C = [   cos(theta3)     sin(theta3)     0;...
                -sin(theta3)    cos(theta3)     0;...
                     0              0           1       ];
         r_body(i,1:3) = r_ACI(i,1:3)*C;
         
    end
end
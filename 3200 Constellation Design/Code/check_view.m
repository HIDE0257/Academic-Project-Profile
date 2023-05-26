%% Check View
function [observable, elevationAngle, cameraAngle,viewMat] = check_view(t,r_body, targets, F, V, Nsc,scNumber,T_benu)
% Function to compute whether a space craft can observe a specified facet
% on the astroid surface

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% r - [x,y,z]' Spacecrafts cartesian position vector in body frame 3x1
% facetNumber - scalar, facet index
% F - matrix of verticies that form each face nx3
% V - matrix of vertex locations in implied body frame
% OUTPUT
% observable - int, 0 for unobservable, 1 for observable
% elevationAngle - scalar, elevation of spacecraft relative to facet plane in radians
% cameraAngle - scalar, angle of facet center relative to camera boresight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Observation requirments:
% local elevation angle must be at least 15 degrees
% center of the facet must be in the feild of view of the camera
% camera always points in the negative radial direction
    X_ACI = zeros(length(t),3) + [-1,0,0];
    X_I = ACI2body(X_ACI,T_benu,t);
    index = 1;
    for j = 1:length(t)
        r = [r_body(j,1) r_body(j,2) r_body(j,3)]';

         for k = 1:length(targets)
            facetNumber = targets(k);
            FacetIndex = F(facetNumber,:); % Verticies corosponding to facet index
            r_center = (sum(V(FacetIndex,:))/3)'; % Center pointiong vetor
            r_view = r-r_center; % Vector from camera to center of facet 
            d = norm(r_view)*1000;
            F_norm = cross(V(F(facetNumber,2),:)-V(F(facetNumber,1),:),V(F(facetNumber,3),:)-V(F(facetNumber,1),:));
            F_norm = F_norm/norm(F_norm);
            H = ceil(dot(X_I(j,:),F_norm));
            elevationAngle = asin(dot(r_view,F_norm)/norm(r_view));
            cameraAngle = acos(dot(r_view,r)/(norm(r_view)*norm(r)));
            FOV = pi/(9*Nsc);
            if elevationAngle >= pi/12 && cameraAngle <= FOV && H == 1
                observable = 1;
                GR = FOV/(2048-(1408/9)*(Nsc-1))*d;
                viewMat(index,1:8) = [facetNumber,t(j),elevationAngle,cameraAngle,FOV,H,d,GR]; 
                index = index+1;
            else
                observable = 0;
            end
           
        end
        
    end
end

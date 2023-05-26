%% Read_OBJ
function [Facets,Vertices] =read_obj(fileName,plot,targets)
% Read in and Plot OBJ File

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name = "Benu"; color = 'y'; dataDelim2 = 'f';
% [facets,Verticies] =read_obj(fileName,Name,color,dataDelim2)
% INPUT
% fileName 
% Name (name of astroid)
% color (color of facets when ploted)
% dataDelim2 (vertical data differentiator)
%
% OUTPUT
%  3D plot of Astroid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Name = "Benu"; 
     color = 'y'; 
     dataDelim2 = 'f';

    FileID = fopen(fileName);
    Contents = textscan(FileID,'%c %f %f %f', 'Delimiter',{''},'CommentStyle',{'#'});
    Type = cell2mat(Contents(1,1));
    Index = find(Type== dataDelim2,1);
    Vector = [cell2mat(Contents(1,2)),cell2mat(Contents(1,3)),cell2mat(Contents(1,4))];
    Vertices = [Vector(1:Index-1,1) Vector(1:Index-1,2) Vector(1:Index-1,3)];
    Facets = [Vector(Index:end,1) Vector(Index:end,2) Vector(Index:end,3)];
    fclose(FileID);
    
    v1 = Vertices(Facets(targets,1),:); 
    v2 = Vertices(Facets(targets,2),:); 
    v3 = Vertices(Facets(targets,3),:); 
    centroid = (v1+v2+v3)./3; 
    if plot == true
        figure
        hold on
        title(Name)
        patch('Faces',Facets,'Vertices',Vertices,'FaceColor',color)
        plot(centroid,'*','r')
        xlabel("X [km]",FontWeight="bold")
        ylabel("Y [km]",FontWeight="bold")
        zlabel("Z [km]",FontWeight="bold")
        view(3)
        hold off
    end
end

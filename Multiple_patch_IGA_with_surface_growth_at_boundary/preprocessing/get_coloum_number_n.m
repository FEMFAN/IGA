function PM = get_coloum_number_n(PM,patches)
% ---------------------------------------------------------------------
% Subroutine get_coloum_number_n.m
% calculation of the coloum for the dofs inside each patch. The PM matrix
% contains the indicies in the gloabl system for the dof inside each patch in a row wise
% manner
%
% Author:           Lennard Langerbein
% Date  :           15.02.2022
%                  
% Input:    PM                                  - global patch matrix
%           patches                             - patch structure
%           
%
% Output:   PM                                  - global patch matrix
%---------------------------------------------------------------------- 

%reset the CM 2 and CM 3 to zero in every time step

%calculation of the coloums. The global system is sorted as the follwoing: At first
%the dofs not inside the interface or next to it for every patch in ascending order KP_n.
%Afterwards the master dofs for every patch starting with dofs next to the
%interface KP_c and lastly the master dofs inside the interface KP_f.
PM(:,2) = 0;
PM(:,3) = 0;
for i = 1:size(PM,1)
    if i ==1
    else
        for j = 1:PM(i,1)-1
            PM(i,2) = PM(i,2) +  length(intersect(patches(j).Int_Dofs(1).KP_n,patches(j).freeDofs));
        end
    end
    PM(i,3) = PM(i,2) + length(intersect(patches(i).Int_Dofs(1).KP_n,patches(i).freeDofs));
end
end
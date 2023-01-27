function IM = get_coloum_number_f(IM,patches,np,continuity_flag)
% ---------------------------------------------------------------------
% Subroutine get_coloum_number_f.m
% calculation of the coloum for the every interface dofs. The IM matrix
% contains the indicies for the gloabl system for the master interface dofs of
% every interface in a row wise manner.
%
% Author:           Lennard Langerbein
% Date  :           15.02.2022
%                  
% Input:    IM                                  - global interface matrix
%           patches                             - patch structure
%           np                                  - number of patches
%           continuity_flag                     - desired continuity
%
% Output:   IM                                  - global interface matrix
%---------------------------------------------------------------------- 

%reset the IM 6 and IM 7 to zero in every time step

%calculation of the coloums. The global system is sorted as the follwoing: At first
%the dofs not inside the interface or next to it for every patch in ascending order KP_n.
%Afterwards the master dofs for every patch starting with dofs next to the
%interface KP_c and lastly the master dofs inside the interface KP_f.
IM(:,6) = 0;
IM(:,7) = 0;
A = 0;
for pl= 1:np
    A= A + length(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs)); 
end
for i = 1:size(IM,1)
    if i ==1
        IM(i,6) = 0;
        if continuity_flag == 1
            IM(i,6) = IM(i,6) + A + length(intersect(patches(IM(i,4)).Int_Dofs(IM(i,5)).KP_f,patches(IM(i,4)).freeDofs));
        elseif continuity_flag == 0
            IM(i,6) = IM(i,6) + A;
        end
        IM(i,7) = IM(i,6) + length(intersect(patches(IM(i,4)).Int_Dofs(IM(i,5)).KP_f,patches(IM(i,4)).freeDofs));
    else
        for j = 1:i-1
            if continuity_flag == 1
                IM(i,6) = IM(i,6) + 2*length(intersect(patches(IM(j,4)).Int_Dofs(IM(j,5)).KP_f,patches(IM(j,4)).freeDofs));
            elseif continuity_flag == 0
                IM(i,6) = IM(i,6) + length(intersect(patches(IM(j,4)).Int_Dofs(IM(j,5)).KP_f,patches(IM(j,4)).freeDofs));
            end
            
        end
        IM(i,6) = IM(i,6) + A+ length(intersect(patches(IM(i,4)).Int_Dofs(IM(i,5)).KP_f,patches(IM(i,4)).freeDofs));
        IM(i,7) = IM(i,6) + length(intersect(patches(IM(i,4)).Int_Dofs(IM(i,5)).KP_f,patches(IM(i,4)).freeDofs));
    end
    
end
end

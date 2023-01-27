function CM = get_coloum_number_c(CM,IM)
% ---------------------------------------------------------------------
% Subroutine get_coloum_number_c.m
% calculation of the coloum for the dofs next to the interfaces. The CM matrix
% contains the indicies in the gloabl system for the master next to the interface dofs  of
% every interface in a row wise manner.
%
% Author:           Lennard Langerbein
% Date  :           15.02.2022
%                  
% Input:    CM                                  - global next to the interface matrix
%           IM                                  - global interface matrix
%           
%
% Output:   CM                                  - global next to the interface matrix
%---------------------------------------------------------------------- 

%reset the CM 2 and CM 3 to zero in every time step

%calculation of the coloums. The global system is sorted as the follwoing: At first
%the dofs not inside the interface or next to it for every patch in ascending order KP_n.
%Afterwards the master dofs for every patch starting with dofs next to the
%interface KP_c and lastly the master dofs inside the interface KP_f.
CM(:,2) = 0;
CM(:,3) = 0;

for i = 1:size(CM,1)
    CM(i,2) = IM(i,6) - IM(i,7) + IM(i,6);
    CM(i,3) = IM(i,7) - IM(i,7) + IM(i,6);
end
function [V] = continuity(u1,u2,KP1,KP2)
% ---------------------------------------------------------------------
% Subroutine continuity.m
% Check for the continuity over the interface. The matrix V contains the
% strains at the interface in a row wise manner. The continuity values are computed over differnce quotient over the interface.
% It only works if both patches are refined equally
%
% Author:           Lennard Langerbein
% Date  :           15.02.2022
%                  
% Input:    u1                                  - displacements of patch 1
%           u2                                  - displacements of patch 2
%           KP1                                 - Control points of patch 1
%           KP1                                 - Control points of patch 2
%           
%
% Output:   V                                   - continuity matrix
%---------------------------------------------------------------------
%initialize all matrices
V = [];
KPc1 = [];
KPf1 = [];
KPc2 = [];
KPf2 = [];
%finde the interface and the corresponding dofs
A = intersect(KP1,KP2,'rows');
for l = 1:size(A,1)
    for i = 1:size(KP1,1)
        if KP1(i,1) == A(l,1) && KP1(i,2) == A(l,2)
            KPf1(end+1) = i*2-1;
            KPf1(end+1) = i*2;
            KPc1(end+1) = i*2-3;
            KPc1(end+1) = i*2-2;
        end
    end
end
for l = 1:size(A,1)
    for i = 1:size(KP2,1)
        if KP2(i,1) == A(l,1) && KP2(i,2) == A(l,2)
            KPf2(end+1) = i*2-1;
            KPf2(end+1) = i*2;
            KPc2(end+1) = i*2+1;
            KPc2(end+1) = i*2+2;
        end
    end
end
%calculation of the continuity matrix as 
for i = 1:2:length(KPc1)
    V(end+1,1) = i;
    V(end,2) = (u1(KPf1(i))-u1(KPc1(i)))/(KP1(round(KPf1(i)/2),1)-(KP1(round(KPc1(i)/2),1)));
    V(end,3) = (u2(KPc2(i))-u2(KPf2(i)))/(KP2(round(KPc2(i)/2),1)-(KP2(round(KPf2(i)/2),1)));
    V(end,4) = abs(V(end,2)-V(end,3));
    %y-value
    V(end,5) = KP1(round(KPf1(i)/2),2);
end

end
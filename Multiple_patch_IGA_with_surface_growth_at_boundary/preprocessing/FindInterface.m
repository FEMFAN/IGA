function [KPf1,KPc1,KPn1,KPf2,KPc2,KPn2] = FindInterface(KP1,KP2,Int,n1,m1,n2,m2,continuity_flag)
% ---------------------------------------------------------------------
% Subroutine FindInterface.m
% find the dofs which belong to the interface, the ones next to the interface 
% and the ones not inside the interface
%
% Author:           Lennard Langerbein
% Date  :           15.02.2022
%                  
% Input:    KP1                                 - control points master patch
%           KP2                                 - control points slave patch
%           Int                                 - interface matrix
%           n1                                  - number of basis function in xi direction of the master patch
%           m1                                  - number of basis function in eta direction of the master patch
%           n2                                  - number of basis function in xi direction of the slave patch
%           n2                                  - number of basis function in eta direction of the slave patch
%           continuity_flag                     - desired continuity in the displacement field
%
% Output:   KPf1                                - dofs inside the master interface
%           KPc1                                - dofs next to the master interface
%           KPn1                                - dofs not in the master interface or next to it 
%           KPf2                                - dofs inside the slave interface
%           KPc2                                - dofs next to the slave interface
%           KPn2                                - dofs not in the slave interface or next to it 
%---------------------------------------------------------------------- 
%initilaize all variables
flag = 0;
KPf1 = [];
KPc1 = [];
KPn1 = [];
KPf2 = [];
KPc2 = [];
KPn2 = [];
%case 1: interface is at the end of the basis function n 
if Int(4) == 0
    %check if the position of the master patch is on the right/upper or on the left/lower part of the interface and switch the number of
    %the basis functions if master patch is on the right/upper side of the
    %interface
    if abs(KP1(n1,1) - KP2(1,1))>1e-8 
        %save the number of basis function
        A=n1;
        B=m1;
        n1 = n2;
        m1 = m2;
        n2=A;
        m2=B;
        flag = 1;
    else
    end
    %sort the dofs
    for i = 0:m1-1
        Indicies =(i*n1+1):((i+1)*n1);
        for j = 1:length(Indicies)
            if j == n1
                KPf1(end+1) = Indicies(j)*2-1;
                KPf1(end+1) = Indicies(j)*2;
            elseif j==n1-1
                if continuity_flag ==1
                    KPc1(end+1) = Indicies(j)*2-1;
                    KPc1(end+1) = Indicies(j)*2;
                else
                    KPn1(end+1) = Indicies(j)*2-1;
                    KPn1(end+1) = Indicies(j)*2;
                    %for the continuity study
                    KPc1(end+1) = Indicies(j)*2-1;
                    KPc1(end+1) = Indicies(j)*2;
                end
            else
                KPn1(end+1) = Indicies(j)*2-1;
                KPn1(end+1) = Indicies(j)*2;
            end
        end
    end
    for i = 0:m2-1
        Indicies =(i*n2+1):((i+1)*n2);
        for j = 1:length(Indicies)
            if j == 1
                KPf2(end+1) = Indicies(j)*2-1;
                KPf2(end+1) = Indicies(j)*2;
            elseif j==2
                if continuity_flag == 1
                    KPc2(end+1) = Indicies(j)*2-1;
                    KPc2(end+1) = Indicies(j)*2;
                else
                    KPn2(end+1) = Indicies(j)*2-1;
                    KPn2(end+1) = Indicies(j)*2;
                    %for the continuity study
                    KPc2(end+1) = Indicies(j)*2-1;
                    KPc2(end+1) = Indicies(j)*2;
                end
            else
                KPn2(end+1) = Indicies(j)*2-1;
                KPn2(end+1) = Indicies(j)*2;
            end
        end
    end
    %change the sorting if flag ==1
    if flag == 1
        A=KPn1;
        B=KPc1;
        C = KPf1;
        KPn1 = KPn2;
        KPc1 = KPc2;
        KPf1 = KPf2;
        KPn2 = A;
        KPc2 = B;
        KPf2 = C;
    else
    end
%case 2: interface is at the end of the basis function m    
elseif Int(4) == 1
    if  abs(KP1(n1*m1-n1+1,2)- KP2(1,2))>1e-8
        %save the number of basis function
        A=n1;
        B=m1;
        n1 = n2;
        m1 = m2;
        n2=A;
        m2=B;
        flag = 1;
    else
    end
    for i = 0:m1-1
        Indicies = (i*n1+1):((i+1)*n1);
        for l = 1:length(Indicies)
            if i == m1-1
                KPf1(end+1) = Indicies(l)*2-1;
                KPf1(end+1) = Indicies(l)*2;
            elseif i == m1-2
                if continuity_flag ==1
                    KPc1(end+1) = Indicies(l)*2-1;
                    KPc1(end+1) = Indicies(l)*2;
                else
                    KPn1(end+1) = Indicies(l)*2-1;
                    KPn1(end+1) = Indicies(l)*2;
                    %for the continuity study
                    KPc1(end+1) = Indicies(l)*2-1;
                    KPc1(end+1) = Indicies(l)*2;
                end
            else
                KPn1(end+1) = Indicies(l)*2-1;
                KPn1(end+1) = Indicies(l)*2;
            end
        end
    end
    for i = 0:m2-1
        Indicies = (i*n2+1):((i+1)*n2);
        for l = 1:length(Indicies)
            if i == 0
                KPf2(end+1) = Indicies(l)*2-1;
                KPf2(end+1) = Indicies(l)*2;
            elseif i == 1
                if continuity_flag ==1
                    KPc2(end+1) = Indicies(l)*2-1;
                    KPc2(end+1) = Indicies(l)*2;
                else
                    KPn2(end+1) = Indicies(l)*2-1;
                    KPn2(end+1) = Indicies(l)*2;
                    %for the continuity study
                    KPc2(end+1) = Indicies(l)*2-1;
                    KPc2(end+1) = Indicies(l)*2;
                end
            else
                KPn2(end+1) = Indicies(l)*2-1;
                KPn2(end+1) = Indicies(l)*2;
            end
        end
    end
    if flag == 1
        A=KPn1;
        B=KPc1;
        C = KPf1;
        KPn1 = KPn2;
        KPc1 = KPc2;
        KPf1 = KPf2;
        KPn2 = A;
        KPc2 = B;
        KPf2 = C;
    else
    end
end
end

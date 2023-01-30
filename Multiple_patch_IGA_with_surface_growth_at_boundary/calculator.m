function [patches,rsn,iter] = calculator(patches,np,iter,tolerance,maxiter,continuity_flag,IM,PM,CM,D1_flag,material_flag,lr_refinement_flag)
% ---------------------------------------------------------------------
% Subroutine Taschenrechner.m
% calculation of the displacements increments in every time step. The used
% method depends on the desired continuity and the refinement
%
% Author:           Lennard Langerbein
% Date  :           15.02.2022
%                  
% Input:    patches                                 - patch strcutur witout the displacements for the current time step
%           np                                      - number of patches
%           iter                                    - current newton iteration
%           tolerance                               - tolerance of the newton loop
%           maxiter                                 - number of maximum iteration
%           continuity_flag                         - desired continuity (0 or 1)
%           IM                                      - global interface matrix
%           CM                                      - global next to the interface matrix
%           PM                                      - global inside the patch matrix
%           D1_flag                                 - flag for 1-dimensional examples
%           lr_refinement_flag                      - flag for local refinement (determines the calculation method in the C^1 case)
%           
%
% Output:   patches                                 - patch structure with the displacements for the current time step
%           rsn                                     - rsn of the newton loop
%           iter                                    - iteration number of the newton loop
%---------------------------------------------------------------------- 

%determine if local refinement was used
local_refinement_flag = any(lr_refinement_flag);

%multiple patch C1 continuity and local refinement
if continuity_flag == 1 && local_refinement_flag == 1
    %check if 2 material example is loaded
    if material_flag== 1 %2 material example -> C^1 continuity is ensured with the colinear method
        %check if the patches are symmetric over the interface precondition
        %for the colinear methods
        for i =1:size(patches(1).Extract_op(1).Csm,1)
           if abs(patches(1).Extract_op(1).Csm(1,1) - patches(1).Extract_op(1).Csm(i,i)) >1e-8
               error('Check symmetry over the interface')
           end
        end
        %begin residuum for all not combined part of the body
        rsd = [];
        for pl = 1:np
            rsd_p = patches(pl).fint(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs)) - patches(pl).fext(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs));
            rsd = [rsd; rsd_p];
        end
        %now all combined parts of the residuum
        rsd_p = patches(1).fint(intersect(patches(1).Int_Dofs(1).KP_c,patches(1).freeDofs)) - patches(1).fext(intersect(patches(1).Int_Dofs(1).KP_c,patches(1).freeDofs)) +...
            patches(2).fint(intersect(patches(2).Int_Dofs(1).KP_c,patches(2).freeDofs)) - patches(2).fext(intersect(patches(2).Int_Dofs(1).KP_c,patches(2).freeDofs));
        rsd = [rsd;rsd_p];
        rsd_p = patches(1).fint(intersect(patches(1).Int_Dofs(1).KP_f,patches(1).freeDofs)) - patches(1).fext(intersect(patches(1).Int_Dofs(1).KP_f,patches(1).freeDofs)) +...
            patches(2).fint(intersect(patches(2).Int_Dofs(1).KP_f,patches(2).freeDofs)) - patches(2).fext(intersect(patches(2).Int_Dofs(1).KP_f,patches(2).freeDofs));
        rsd = [rsd;rsd_p];
        A = size(rsd,1);
        %determine the norm of the residuum
        rsn = norm(rsd);
        fprintf(1, ' %4d, residuum_norm= %e\n', iter, rsn);
        if rsn > tolerance
                %beginn with the global system matrix
                J = zeros(A);
                %determine help varables for the sake of better reading
                N1 = length(intersect(patches(1).Int_Dofs(1).KP_n,patches(1).freeDofs)); %number of free Dofs inside patch 1
                N2 = length(intersect(patches(2).Int_Dofs(1).KP_n,patches(2).freeDofs)); %number of free Dofs inside patch 2
                C1 = length(intersect(patches(1).Int_Dofs(1).KP_c,patches(1).freeDofs)); %number of free Dofs near the interface of patch 1
                C2 = length(intersect(patches(2).Int_Dofs(1).KP_c,patches(2).freeDofs)); %number of free Dofs near the interface of patch 2
                %determine the factor C
                c = patches(1).Extract_op(1).Csm/(patches(1).Extract_op(1).Csm+eye(size(patches(1).Extract_op(1).Csm,1)));
                %begin with the not combined parts of the patches
                %first row (patch 1)
                J(1:N1,1:N1)                    = patches(1).K_m(1).K_nn;
                J(1:N1,N1+N2+1:N1+N2+C1)        = patches(1).K_m(1).K_nc +patches(1).K_m(1).K_nf*c;
                J(1:N1,N1+N2+C1+1:N1+N2+C1+C2)  = patches(1).K_m(1).K_nf*c;
                %second row (patch 2)
                J(N1+1:N1+N2,N1+1:N1+N2)              = patches(2).K_m(1).K_nn;
                J(N1+1:N1+N2,N1+N2+1:N1+N2+C1)        = patches(2).K_m(1).K_nf*c;
                J(N1+1:N1+N2,N1+N2+C1+1:N1+N2+C1+C2)  = patches(2).K_m(1).K_nc + patches(2).K_m(1).K_nf*c;
                %next the combined parts of the system, beginning with the
                %close to the interface part
                J(N1+N2+1:N1+N2+C1,1:N1)                    = patches(1).K_m(1).K_cn;
                J(N1+N2+1:N1+N2+C1,N1+1:N1+N2)              = patches(2).K_m(1).K_cn;
                J(N1+N2+1:N1+N2+C1,N1+N2+1:N1+N2+C1)        = patches(1).K_m(1).K_cc + (patches(1).K_m(1).K_cf+patches(2).K_m(1).K_cf)*c;
                J(N1+N2+1:N1+N2+C1,N1+N2+C1+1:N1+N2+C1+C2)  = patches(2).K_m(1).K_cc + (patches(1).K_m(1).K_cf+patches(2).K_m(1).K_cf)*c;
                %the interface part
                J(N1+N2+C1+1:N1+N2+C1+C2,1:N1)                    = patches(1).K_m(1).K_fn;
                J(N1+N2+C1+1:N1+N2+C1+C2,N1+1:N1+N2)              = patches(2).K_m(1).K_fn;
                J(N1+N2+C1+1:N1+N2+C1+C2,N1+N2+1:N1+N2+C1)        = patches(1).K_m(1).K_fc + (patches(1).K_m(1).K_ff+patches(2).K_m(1).K_ff)*c;
                J(N1+N2+C1+1:N1+N2+C1+C2,N1+N2+C1+1:N1+N2+C1+C2)  = patches(2).K_m(1).K_fc + (patches(1).K_m(1).K_ff+patches(2).K_m(1).K_ff)*c;
                %determine the displacements increments
                du = -J\rsd;
                %write the dislacements inside the patch structure
                patches(1).u(intersect(patches(1).Int_Dofs(1).KP_n,patches(1).freeDofs)) = patches(1).u(intersect(patches(1).Int_Dofs(1).KP_n,patches(1).freeDofs)) + du(1:N1);
                patches(2).u(intersect(patches(2).Int_Dofs(1).KP_n,patches(2).freeDofs)) = patches(2).u(intersect(patches(2).Int_Dofs(1).KP_n,patches(2).freeDofs)) + du(N1+1:N1+N2);
                patches(1).u(intersect(patches(1).Int_Dofs(1).KP_c,patches(1).freeDofs)) = patches(1).u(intersect(patches(1).Int_Dofs(1).KP_c,patches(1).freeDofs)) + du(N1+N2+1:N1+N2+C1);
                patches(2).u(intersect(patches(2).Int_Dofs(1).KP_c,patches(2).freeDofs)) = patches(2).u(intersect(patches(2).Int_Dofs(1).KP_c,patches(2).freeDofs)) + du(N1+N2+C1+1:N1+N2+C1+C2);
                
                patches(1).u(intersect(patches(1).Int_Dofs(1).KP_f,patches(1).freeDofs)) = c* (patches(1).u(intersect(patches(1).Int_Dofs(1).KP_c,patches(1).freeDofs)) + patches(2).u(intersect(patches(2).Int_Dofs(1).KP_c,patches(2).freeDofs)) );
                patches(2).u(intersect(patches(2).Int_Dofs(1).KP_f,patches(2).freeDofs)) = patches(1).u(intersect(patches(1).Int_Dofs(1).KP_f,patches(1).freeDofs));
                
                iter = iter +1;
                if iter == maxiter
                    close all
                    error('No convergence in Newton scheme!')
                end
           end
    else %1 material examples -> C^1 continuity is ensured with the VUKIMS method
        %begin residuum
        %initialise the residuum for all not combined parts of the body
        rsd = [];
        for pl = 1:np
            rsd_p = patches(pl).fint(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs)) - patches(pl).fext(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs));
            rsd = [rsd; rsd_p];
        end
        A = size(rsd,1);
        %now all combined parts of the residuum for every interface
        for pl = 1:size(IM,1)
            ma=IM(pl,4);%master patch
            for i = 1:size(patches(ma).Int,1)
                if pl == patches(ma).Int(i,1)
                    s =  patches(ma).Int(i,5); %slave patch number
                    %position inside the slave patch interface matrix
                    for l = 1:size(patches(s).Int,1)
                        if pl == patches(s).Int(l,1)
                            p = l; %position inside the slave interface matrix
                        else
                        end
                    end
                    %calculate the extension operators for the current Interface
                    if D1_flag == 1
                        Csm = patches(s).Extract_op(p).Csm_1D;
                        Ts1m1 = pinv(patches(s).Extract_op(p).A1_1D) *patches(ma).Extract_op(i).A1_1D;
                        Ts2m2 = -pinv(patches(s).Extract_op(p).A2_1D) *Csm*patches(ma).Extract_op(i).A2_1D;
                    else
                        Csm = patches(s).Extract_op(p).Csm;
                        Ts1m1 = pinv(patches(s).Extract_op(p).A1) *patches(ma).Extract_op(i).A1;
                        Ts2m2 = -pinv(patches(s).Extract_op(p).A2) *Csm*patches(ma).Extract_op(i).A2;
                    end
                    %first the residuum for the c part
                    rsd_p = patches(ma).fint(intersect(patches(ma).Int_Dofs(i).KP_c,patches(ma).freeDofs)) - patches(ma).fext(intersect(patches(ma).Int_Dofs(i).KP_c,patches(ma).freeDofs)) +...
                        Ts2m2'*patches(s).fint(intersect(patches(s).Int_Dofs(p).KP_c,patches(s).freeDofs)) - Ts2m2'*patches(s).fext(intersect(patches(s).Int_Dofs(p).KP_c,patches(s).freeDofs));
                    rsd = [rsd; rsd_p];
                    %second the residuum for the f part
                    rsd_p = patches(s).fint(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs)) - patches(s).fext(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs));
                    
                    rsd_p = patches(ma).fint(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) - patches(ma).fext(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) +...
                        Ts1m1'*patches(s).fint(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs)) - Ts1m1'*patches(s).fext(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs));
                    rsd = [rsd; rsd_p];
                end
            end
            B = size(rsd,1);
        end
        %end residuum
        rsn = norm(rsd);
        fprintf(1, ' %4d, residuum_norm= %e\n', iter, rsn);
        if rsn > tolerance
            %begin stiffness matrix
            J1 = zeros(A,size(rsd,1));
            J2 = zeros(B-A,size(rsd,1));
            %first fill the uncombined parts of every patch in ascending order
            A = 0; %counting variable
            %first the K_nn matrices of every patch
            for pl = 1:np
                NP = length(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs));%help variable for the sake of better reading
                J1(1+A:NP+A,1+A:NP+A) = patches(pl).K_m(1).K_nn;
                
                %second the K matrices from all n CP to all other interfaces
                %(K_nci and K_nfi)
                for i = 1:size(patches(pl).Int,1)
                    %first the K_nci matrices
                    %check if the current patch is the master patch of the
                    %interface and get the coloums of this matrix, also get the
                    %position inside the masters interface matrix
                    [alpha,omega,ma,pci] = get_coloum_c(patches(pl).Int(i,1),IM,CM);
                    %calculate the expansion operator Ts2m2
                    if ma==pl %the current patch is the master patch of this interface
                        J1(1+A:NP+A,1+alpha:omega) = patches(pl).K_m(i).K_nc;
                        %second the K_nfi matrices
                        %get the coloums and fill in the matrix
                        [alpha,omega,~,~] = get_coloum_f(patches(pl).Int(i,1),IM);
                        J1(1+A:NP+A,1+alpha:omega) = patches(pl).K_m(i).K_nf;
                    elseif ma~=pl %the current patch is not the master patch of this interface
                        if D1_flag == 1
                            Csm = patches(pl).Extract_op(i).Csm_1D;
                            Ts1m1 = pinv(patches(pl).Extract_op(i).A1_1D) *patches(ma).Extract_op(pci).A1_1D;
                            Ts2m1 = pinv(patches(pl).Extract_op(i).A2_1D)*(patches(pl).Extract_op(i).A1_1D*pinv(patches(pl).Extract_op(i).A1_1D)+Csm)*patches(ma).Extract_op(pci).A1_1D;
                            Ts2m2 = -pinv(patches(pl).Extract_op(i).A2_1D) *Csm*patches(ma).Extract_op(pci).A2_1D;
                        else
                            Csm = patches(pl).Extract_op(i).Csm;
                            Ts1m1 = pinv(patches(pl).Extract_op(i).A1) *patches(ma).Extract_op(pci).A1;
                            Ts2m1 = pinv(patches(pl).Extract_op(i).A2)*(patches(pl).Extract_op(i).A1*pinv(patches(pl).Extract_op(i).A1)+Csm)*patches(ma).Extract_op(pci).A1;
                            Ts2m2 = -pinv(patches(pl).Extract_op(i).A2) *Csm*patches(ma).Extract_op(pci).A2;
                        end
                        J1(1+A:NP+A,1+alpha:omega) = patches(pl).K_m(i).K_nc * Ts2m2;
                        %second the K_nfi matrices
                        %get the coloums and fill in the matrix
                        [alpha,omega,~,~] = get_coloum_f(patches(pl).Int(i,1),IM);
                        J1(1+A:NP+A,1+alpha:omega) = patches(pl).K_m(i).K_nc * Ts2m1 + patches(pl).K_m(i).K_nf * Ts1m1;
                    end
                    
                    
                end
                A = A + NP;
            end
            %second all the combined parts going over the interface
            %starting with c part of the interface
            A = 0; %counting variable
            for il = 1:size(IM,1)%loop over all interfaces
                m = IM(il,4); %master patch of the current interface
                for i = 1:size(patches(m).Int,1) %finding the interface inside the patch inteface matrix
                    if il == patches(m).Int(i,1)
                        s = patches(m).Int(i,5); %slave patch
                        NF = length(intersect(patches(m).Int_Dofs(i).KP_c,patches(m).freeDofs));%help variable for the sake of better reading
                        %first fill in the two K_cni matrices
                        %first the master K_cni matrix
                        %find the coloumns of the master patch
                        [alpha,omega] = get_coloum_n(PM,m);
                        J2(1+A:NF+A,1+alpha:omega) = patches(m).K_m(i).K_cn;
                        %second the slave K_cni matrix
                        %find the coloumns of the slave patch
                        [alpha,omega] = get_coloum_n(PM,s);
                        %find the position inside the slaves interface matrix
                        for l = 1:size(patches(s).Int,1)
                            if il == patches(s).Int(l,1)
                                p = l; %position inside the slave interface matrix
                            else
                            end
                        end
                        %calculate the extension operators
                        if D1_flag == 1
                            Csm = patches(s).Extract_op(p).Csm_1D;
                            Ts1m1 = pinv(patches(s).Extract_op(p).A1_1D)*patches(m).Extract_op(i).A1_1D;
                            Ts2m1 = pinv(patches(s).Extract_op(p).A2_1D)*(patches(s).Extract_op(p).A1_1D*pinv(patches(s).Extract_op(p).A1_1D)+Csm)*patches(m).Extract_op(i).A1_1D;
                            Ts2m2 = -pinv(patches(s).Extract_op(p).A2_1D) *Csm*patches(m).Extract_op(i).A2_1D;
                        else
                            Csm = patches(s).Extract_op(p).Csm;
                            Ts1m1 = pinv(patches(s).Extract_op(p).A1)*patches(m).Extract_op(i).A1;
                            Ts2m1 = pinv(patches(s).Extract_op(p).A2)*(patches(s).Extract_op(p).A1*pinv(patches(s).Extract_op(p).A1)+Csm)*patches(m).Extract_op(i).A1;
                            Ts2m2 = -pinv(patches(s).Extract_op(p).A2) *Csm*patches(m).Extract_op(i).A2;
                        end
                        %fill in the matrix
                        J2(1+A:NF+A,1+alpha:omega) = Ts2m2'*patches(s).K_m(p).K_cn;
                        %next the combined K_cc matrices on the current interface
                        %find the indicies first
                        [alpha,omega,~,~] = get_coloum_c(il,IM,CM);
                        %fill in the matrices
                        J2(1+A:NF+A,1+alpha:omega) = patches(m).K_m(i).K_cc +  Ts2m2'*patches(s).K_m(p).K_cc*Ts2m2;
                        %after this the K_cf matrices on the current interface
                        %find the indicies first
                        [alpha,omega,~,~] = get_coloum_f(il,IM);
                        %fill in the matrices
                        J2(1+A:NF+A,1+alpha:omega) = patches(m).K_m(i).K_cf +  Ts2m2'*(patches(s).K_m(p).K_cc*Ts2m1+patches(s).K_m(p).K_cf*Ts1m1);
                        %now all stiffnes matrices between the master c points
                        %and other c and f points of other interfaces
                        for l = 1:size(patches(m).Int,1)
                            if il ~= patches(m).Int(l,1)
                                %check if the current patch is the master of this
                                %patch
                                [alpha,omega,ma,pci] = get_coloum_c(patches(m).Int(l,1),IM,CM);
                                if ma ==m %the master patch is also the master patch of the other interface
                                    %fill in the matrices
                                    %first the K_cci matrix
                                    K_cci = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_c,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_c,patches(m).freeDofs));
                                    %then the K_cfi matrix
                                    K_cfi = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_c,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_f,patches(m).freeDofs));
                                    %fill in the matrices
                                    %first the K_cci matrices
                                    %coloums already calculated
                                    J2(1+A:NF+A,1+alpha:omega) =  K_cci;
                                    %second the K_cfi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(m).Int(l,1),IM);
                                    %fill in the matrix
                                    J2(1+A:NF+A,1+alpha:omega) = K_cfi;
                                    
                                elseif ma~=m %the master patch is the slave patch of the other interface
                                    %first calculate the extension operator
                                    if D1_flag == 1
                                        Csm = patches(m).Extract_op(l).Csm_1D;
                                        Ts1m1 = pinv(patches(m).Extract_op(l).A1_1D)*patches(ma).Extract_op(pci).A1_1D;
                                        Ts2m1 = pinv(patches(m).Extract_op(l).A2_1D)*(patches(m).Extract_op(l).A1_1D*pinv(patches(m).Extract_op(l).A1_1D)+Csm)*patches(ma).Extract_op(pci).A1_1D;
                                        Ts2m2 = -pinv(patches(m).Extract_op(l).A2_1D) *Csm*patches(ma).Extract_op(pci).A2_1D;
                                    else
                                        Csm = patches(m).Extract_op(l).Csm;
                                        Ts1m1 = pinv(patches(m).Extract_op(l).A1)*patches(ma).Extract_op(pci).A1;
                                        Ts2m1 = pinv(patches(m).Extract_op(l).A2)*(patches(m).Extract_op(l).A1*pinv(patches(m).Extract_op(l).A1)+Csm)*patches(ma).Extract_op(pci).A1;
                                        Ts2m2 = -pinv(patches(m).Extract_op(l).A2) *Csm*patches(ma).Extract_op(pci).A2;
                                    end
                                    %fill in the matrices
                                    %first the K_cci matrices
                                    %coloums already calculated
                                    %get the stiffnes matrix
                                    K_cci = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_c,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_c,patches(m).freeDofs));
                                    J2(1+A:NF+A,1+alpha:omega) =  K_cci*Ts2m2;
                                    %second the K_cfi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(m).Int(l,1),IM);
                                    %get the stiffnes matrix
                                    K_cfi = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_c,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_f,patches(m).freeDofs));
                                    %fill in the matrix
                                    J2(1+A:NF+A,1+alpha:omega) = K_cci*Ts2m1 + K_cfi*Ts1m1;
                                    
                                end
                                
                            end
                        end
                        %now all stiffnes matrices between the slave c points
                        %and other c and f points of other interfaces
                        %recalculate the extension operator of this interface
                        if D1_flag == 1
                            Csm = patches(s).Extract_op(p).Csm_1D;
                            Ts2m2I = -pinv(patches(s).Extract_op(p).A2_1D) *Csm*patches(m).Extract_op(i).A2_1D;
                        else
                            Csm = patches(s).Extract_op(p).Csm;
                            Ts2m2I = -pinv(patches(s).Extract_op(p).A2) *Csm*patches(m).Extract_op(i).A2;
                        end
                        for l = 1:size(patches(s).Int,1)
                            if il ~= patches(s).Int(l,1)
                                %check if the current patch is the master of this
                                %interface
                                [alpha,omega,ma,pci] = get_coloum_c(patches(s).Int(l,1),IM,CM);
                                if ma ==s %the slave patch is the master patch of the other interface
                                    %fill in the matrices
                                    %first the K_cci matrix
                                    K_cci = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_c,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_c,patches(s).freeDofs));
                                    %then the K_cfi matrix
                                    K_cfi = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_c,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_f,patches(s).freeDofs));
                                    %fill in the matrices
                                    %first the K_cci matrices
                                    %coloums already calculated
                                    J2(1+A:NF+A,1+alpha:omega) =  Ts2m2I'*K_cci;
                                    %second the K_cfi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(s).Int(l,1),IM);
                                    %fill in the matrix
                                    J2(1+A:NF+A,1+alpha:omega) = Ts2m2I'*K_cfi;
                                    
                                elseif ma~=m %the slave patch is also the slave patch of the other interface
                                    %first calculate the extension operator
                                    if D1_flag == 1
                                        Csm = patches(s).Extract_op(l).Csm_1D;
                                        Ts1m1 = pinv(patches(s).Extract_op(l).A1_1D)*patches(ma).Extract_op(pci).A1_1D;
                                        Ts2m1 = pinv(patches(s).Extract_op(l).A2_1D)*(patches(s).Extract_op(l).A1_1D*pinv(patches(s).Extract_op(l).A1_1D)+Csm)*patches(ma).Extract_op(pci).A1_1D;
                                        Ts2m2 = -pinv(patches(s).Extract_op(l).A2_1D) *Csm*patches(ma).Extract_op(pci).A2_1D;
                                    else
                                        Csm = patches(s).Extract_op(l).Csm;
                                        Ts1m1 = pinv(patches(s).Extract_op(l).A1)*patches(ma).Extract_op(pci).A1;
                                        Ts2m1 = pinv(patches(s).Extract_op(l).A2)*(patches(s).Extract_op(l).A1*pinv(patches(s).Extract_op(l).A1)+Csm)*patches(ma).Extract_op(pci).A1;
                                        Ts2m2 = -pinv(patches(s).Extract_op(l).A2) *Csm*patches(ma).Extract_op(pci).A2;
                                    end
                                    %fill in the matrices
                                    %first the K_cci matrices
                                    %coloums already calculated
                                    %get the stiffnes matrix
                                    K_cci = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_c,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_c,patches(s).freeDofs));
                                    J2(1+A:NF+A,1+alpha:omega) =   Ts2m2I'*K_cci*Ts2m2;
                                    %second the K_cfi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(s).Int(l,1),IM);
                                    %get the stiffnes matrix
                                    K_cfi = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_c,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_f,patches(s).freeDofs));
                                    %fill in the matrix
                                    J2(1+A:NF+A,1+alpha:omega) =  Ts2m2I'*(K_cci*Ts2m1 + K_cfi*Ts1m1);
                                    
                                end
                            end
                            
                        end
                        %now the f part of the interface
                        %first fill in the two K_fni matrices
                        %first the master K_fni matrix
                        %find the coloumns of the master patch
                        [alpha,omega] = get_coloum_n(PM,m);
                        J2(1+A+NF:2*NF+A,1+alpha:omega) = patches(m).K_m(i).K_fn;
                        %second the slave K_fni matrix
                        %find the coloumns of the slave patch
                        [alpha,omega] = get_coloum_n(PM,s);
                        %find the position inside the slaves interface matrix
                        for l = 1:size(patches(s).Int,1)
                            if il == patches(s).Int(l,1)
                                p = l; %position inside the slave interface matrix
                            else
                            end
                        end
                        %calculate the extension operators
                        if D1_flag == 1
                            Csm = patches(s).Extract_op(p).Csm_1D;
                            Ts1m1 = pinv(patches(s).Extract_op(p).A1_1D)*patches(m).Extract_op(i).A1_1D;
                            Ts2m1 = pinv(patches(s).Extract_op(p).A2_1D)*(patches(s).Extract_op(p).A1_1D*pinv(patches(s).Extract_op(p).A1_1D)+Csm)*patches(m).Extract_op(i).A1_1D;
                            Ts2m2 = -pinv(patches(s).Extract_op(p).A2_1D) *Csm*patches(m).Extract_op(i).A2_1D;
                        else
                            Csm = patches(s).Extract_op(p).Csm;
                            Ts1m1 = pinv(patches(s).Extract_op(p).A1)*patches(m).Extract_op(i).A1;
                            Ts2m1 = pinv(patches(s).Extract_op(p).A2)*(patches(s).Extract_op(p).A1*pinv(patches(s).Extract_op(p).A1)+Csm)*patches(m).Extract_op(i).A1;
                            Ts2m2 = -pinv(patches(s).Extract_op(p).A2) *Csm*patches(m).Extract_op(i).A2;
                        end
                        %fill in the matrix
                        J2(1+A+NF:2*NF+A,1+alpha:omega) = Ts1m1'*patches(s).K_m(p).K_fn;
                        %next the combined K_ff matrices on the current interface
                        %find the indicies first
                        [alpha,omega,~,~] = get_coloum_f(il,IM);
                        %fill in the matrices
                        J2(1+A+NF:2*NF+A,1+alpha:omega) = patches(m).K_m(i).K_ff +  Ts1m1'*(patches(s).K_m(p).K_fc*Ts2m1 + patches(s).K_m(p).K_ff*Ts1m1);
                        %after this the K_fc matrices on the current interface
                        %find the indicies first
                        [alpha,omega,~,~] = get_coloum_c(il,IM,CM);
                        %fill in the matrices
                        J2(1+A+NF:2*NF+A,1+alpha:omega) = patches(m).K_m(i).K_fc +  Ts1m1'*patches(s).K_m(p).K_fc*Ts2m2;
                        %now all stiffnes matrices between the master f points
                        %and other c and f points of other interfaces
                        for l = 1:size(patches(m).Int,1)
                            if il ~= patches(m).Int(l,1)
                                %check if the current patch is the master of this
                                %patch
                                [alpha,omega,ma,pci] = get_coloum_c(patches(m).Int(l,1),IM,CM);
                                if ma ==m %the master patch is also the master patch of the other interface
                                    %fill in the matrices
                                    %first the K_fci matrix
                                    K_fci = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_f,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_c,patches(m).freeDofs));
                                    %then the K_ffi matrix
                                    K_ffi = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_f,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_f,patches(m).freeDofs));
                                    %fill in the matrices
                                    %first the K_fci matrices
                                    %coloums already calculated
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) =  K_fci;
                                    %second the K_ffi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(m).Int(l,1),IM);
                                    %fill in the matrix
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) = K_ffi;
                                    
                                elseif ma~=m %the master patch is the slave patch of the other interface
                                    %first calculate the extension operator
                                    if D1_flag == 1
                                        Csm = patches(m).Extract_op(l).Csm_1D;
                                        Ts1m1 = pinv(patches(m).Extract_op(l).A1_1D)*patches(ma).Extract_op(pci).A1_1D;
                                        Ts2m1 = pinv(patches(m).Extract_op(l).A2_1D)*(patches(m).Extract_op(l).A1_1D*pinv(patches(m).Extract_op(l).A1_1D)+Csm)*patches(ma).Extract_op(pci).A1_1D;
                                        Ts2m2 = -pinv(patches(m).Extract_op(l).A2_1D) *Csm*patches(ma).Extract_op(pci).A2_1D;
                                    else
                                        Csm = patches(m).Extract_op(l).Csm;
                                        Ts1m1 = pinv(patches(m).Extract_op(l).A1)*patches(ma).Extract_op(pci).A1;
                                        Ts2m1 = pinv(patches(m).Extract_op(l).A2)*(patches(m).Extract_op(l).A1*pinv(patches(m).Extract_op(l).A1)+Csm)*patches(ma).Extract_op(pci).A1;
                                        Ts2m2 = -pinv(patches(m).Extract_op(l).A2) *Csm*patches(ma).Extract_op(pci).A2;
                                    end
                                    %fill in the matrices
                                    %first the K_cci matrices
                                    %coloums already calculated
                                    %get the stiffnes matrix
                                    K_fci = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_f,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_c,patches(m).freeDofs));
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) =  K_fci*Ts2m2;
                                    %second the K_ffi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(m).Int(l,1),IM);
                                    %get the stiffnes matrix
                                    K_ffi = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_f,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_f,patches(m).freeDofs));
                                    %fill in the matrix
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) = K_fci*Ts2m1 + K_ffi*Ts1m1;
                                    
                                end
                                
                            end
                        end
                        %now all stiffnes matrices between the slave f points
                        %and other c and f points of other interfaces
                        %recalculate the extension operator of this interface
                        if D1_flag == 1
                            Csm = patches(s).Extract_op(p).Csm_1D;
                            Ts2m2I = -pinv(patches(s).Extract_op(p).A2_1D) *Csm*patches(m).Extract_op(i).A2_1D;
                        else
                            Csm = patches(s).Extract_op(p).Csm;
                            Ts2m2I = -pinv(patches(s).Extract_op(p).A2) *Csm*patches(m).Extract_op(i).A2;
                        end
                        for l = 1:size(patches(s).Int,1)
                            if il ~= patches(s).Int(l,1)
                                %check if the current patch is the master of this
                                %patch
                                [alpha,omega,ma,pci] = get_coloum_c(patches(s).Int(l,1),IM,CM);
                                if ma ==s %the slave patch is the master patch of the other interface
                                    %fill in the matrices
                                    %first the K_fci matrix
                                    K_fci = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_c,patches(s).freeDofs));
                                    %then the K_ffi matrix
                                    K_ffi = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_f,patches(s).freeDofs));
                                    %fill in the matrices
                                    %first the K_fci matrices
                                    %coloums already calculated
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) =  Ts2m2I'*K_fci;
                                    %second the K_ffi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(s).Int(l,1),IM);
                                    %fill in the matrix
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) = Ts2m2I'*K_ffi;
                                    
                                elseif ma~=s %the slave patch is also the slave patch of the other interface
                                    %first calculate the extension operators
                                    if D1_flag == 1
                                        Csm = patches(s).Extract_op(l).Csm_1D;
                                        Ts1m1 = pinv(patches(s).Extract_op(l).A1_1D)*patches(ma).Extract_op(pci).A1_1D;
                                        Ts2m1 = pinv(patches(s).Extract_op(l).A2_1D)*(patches(s).Extract_op(l).A1_1D*pinv(patches(s).Extract_op(l).A1_1D)+Csm)*patches(ma).Extract_op(pci).A1_1D;
                                        Ts2m2 = -pinv(patches(s).Extract_op(l).A2_1D) *Csm*patches(ma).Extract_op(pci).A2_1D;
                                    else
                                        Csm = patches(s).Extract_op(l).Csm;
                                        Ts1m1 = pinv(patches(s).Extract_op(l).A1)*patches(ma).Extract_op(pci).A1;
                                        Ts2m1 = pinv(patches(s).Extract_op(l).A2)*(patches(s).Extract_op(l).A1*pinv(patches(s).Extract_op(l).A1)+Csm)*patches(ma).Extract_op(pci).A1;
                                        Ts2m2 = -pinv(patches(s).Extract_op(l).A2) *Csm*patches(ma).Extract_op(pci).A2;
                                    end
                                    %fill in the matrices
                                    %first the K_fci matrices
                                    %coloums already calculated
                                    %get the stiffnes matrix
                                    K_fci = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_c,patches(s).freeDofs));
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) =   Ts2m2I'*K_fci*Ts2m2;
                                    %second the K_ffi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(s).Int(l,1),IM);
                                    %get the stiffnes matrix
                                    K_ffi = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_f,patches(s).freeDofs));
                                    %fill in the matrix
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) =  Ts2m2I'*(K_fci*Ts2m1 + K_ffi*Ts1m1);
                                end
                            end
                        end
                        
                    end
                end
                A = A+ 2*NF;
            end
            J = [J1; J2];
            %calculate the displacement increments
            du = -J\rsd;
            %write the displacments into the patch structure
            %first the u_n displacements for every patch
            for pl = 1:np
                %get the indicies of the patch
                [alpha,omega] = get_coloum_n(PM,pl);
                patches(pl).u(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs)) = patches(pl).u(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs)) + du(alpha+1:omega);
            end
            for pl = 1:size(IM,1)
                ma = IM(pl,4);
                for i = 1:size(patches(ma).Int,1)
                    if pl==patches(ma).Int(i,1)
                        %get the indicies of the master patch displacements
                        [alpha,omega,m,pci] = get_coloum_f(patches(ma).Int(i,1),IM);
                        %write the master displacements inside the patch
                        %structure
                        patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) = ...
                            patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) + du(alpha+1:omega);
                        %get the indicies of the master patch displacements
                        [alpha,omega,m,pci] = get_coloum_c(patches(ma).Int(i,1),IM,CM);
                        %write the master displacements inside the patch
                        %structure
                        patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_c,patches(ma).freeDofs)) = ...
                            patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_c,patches(ma).freeDofs)) + du(alpha+1:omega);
                        %get the slave patch and the indice inside the slave
                        %interface matrix
                        s = patches(ma).Int(i,5);
                        for l = 1:size(patches(s).Int,1)
                            if pl == patches(s).Int(l,1)
                                p =l;
                            end
                        end
                        %calculate the relation operator T and the slave
                        %displacements inside the interface
                        if D1_flag == 1
                            Csm = patches(s).Extract_op(p).Csm_1D;
                            Ts1m1 = pinv(patches(s).Extract_op(p).A1_1D)*patches(ma).Extract_op(i).A1_1D;
                            Ts2m1 = pinv(patches(s).Extract_op(p).A2_1D)*(patches(s).Extract_op(p).A1_1D*pinv(patches(s).Extract_op(p).A1_1D)+Csm)*patches(ma).Extract_op(i).A1_1D;
                            Ts2m2 = -pinv(patches(s).Extract_op(p).A2_1D) *Csm*patches(ma).Extract_op(i).A2_1D;
                        else
                            Csm = patches(s).Extract_op(p).Csm;
                            Ts1m1 = pinv(patches(s).Extract_op(p).A1)*patches(ma).Extract_op(i).A1;
                            Ts2m1 = pinv(patches(s).Extract_op(p).A2)*(patches(s).Extract_op(p).A1*pinv(patches(s).Extract_op(p).A1)+Csm)*patches(ma).Extract_op(i).A1;
                            Ts2m2 = -pinv(patches(s).Extract_op(p).A2) *Csm*patches(ma).Extract_op(i).A2;
                        end
                        patches(s).u(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs)) = ...
                            Ts1m1 *  patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs));
                        patches(s).u(intersect(patches(s).Int_Dofs(p).KP_c,patches(s).freeDofs)) = ...
                            Ts2m1 *  patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) +...
                            Ts2m2 * patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_c,patches(ma).freeDofs));
                    end
                end
            end
            iter = iter +1;
            if iter == maxiter
                close all
                error('No convergence in Newton scheme!')
            end
        end
    end
end

%multiple patch C1 continuity and no local refinement
if continuity_flag == 1 && local_refinement_flag == 0
    %check if 2 material example is loaded
    if material_flag== 1 %2 material example -> C^1 continuity is ensured with the colinear method
        %check if the patches are symmetric over the interface precondition
        %for the colinear methods
        for i =1:size(patches(1).Extract_op(1).Csm,1)
           if abs(patches(1).Extract_op(1).Csm(1,1) - patches(1).Extract_op(1).Csm(i,i)) >1e-8
               error('Check symmetry over the interface')
           end
        end
        %begin residuum for all not combined part of the body
        rsd = [];
        for pl = 1:np
            rsd_p = patches(pl).fint(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs)) - patches(pl).fext(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs));
            rsd = [rsd; rsd_p];
        end
        %now all combined parts of the residuum
        rsd_p = patches(1).fint(intersect(patches(1).Int_Dofs(1).KP_c,patches(1).freeDofs)) - patches(1).fext(intersect(patches(1).Int_Dofs(1).KP_c,patches(1).freeDofs)) +...
            patches(2).fint(intersect(patches(2).Int_Dofs(1).KP_c,patches(2).freeDofs)) - patches(2).fext(intersect(patches(2).Int_Dofs(1).KP_c,patches(2).freeDofs));
        rsd = [rsd;rsd_p];
        rsd_p = patches(1).fint(intersect(patches(1).Int_Dofs(1).KP_f,patches(1).freeDofs)) - patches(1).fext(intersect(patches(1).Int_Dofs(1).KP_f,patches(1).freeDofs)) +...
            patches(2).fint(intersect(patches(2).Int_Dofs(1).KP_f,patches(2).freeDofs)) - patches(2).fext(intersect(patches(2).Int_Dofs(1).KP_f,patches(2).freeDofs));
        rsd = [rsd;rsd_p];
        A = size(rsd,1);
        %determine the norm of the residuum
        rsn = norm(rsd);
        fprintf(1, ' %4d, residuum_norm= %e\n', iter, rsn);
        if rsn > tolerance
                %beginn with the global system matrix
                J = zeros(A);
                %determine help varables for the sake of better reading
                N1 = length(intersect(patches(1).Int_Dofs(1).KP_n,patches(1).freeDofs)); %number of free Dofs inside patch 1
                N2 = length(intersect(patches(2).Int_Dofs(1).KP_n,patches(2).freeDofs)); %number of free Dofs inside patch 2
                C1 = length(intersect(patches(1).Int_Dofs(1).KP_c,patches(1).freeDofs)); %number of free Dofs near the interface of patch 1
                C2 = length(intersect(patches(2).Int_Dofs(1).KP_c,patches(2).freeDofs)); %number of free Dofs near the interface of patch 2
                %determine the factor C
                c = patches(1).Extract_op(1).Csm/(patches(1).Extract_op(1).Csm+eye(size(patches(1).Extract_op(1).Csm,1)));
                %begin with the not combined parts of the patches
                %first row (patch 1)
                J(1:N1,1:N1)                    = patches(1).K_m(1).K_nn;
                J(1:N1,N1+N2+1:N1+N2+C1)        = patches(1).K_m(1).K_nc +patches(1).K_m(1).K_nf*c;
                J(1:N1,N1+N2+C1+1:N1+N2+C1+C2)  = patches(1).K_m(1).K_nf*c;
                %second row (patch 2)
                J(N1+1:N1+N2,N1+1:N1+N2)              = patches(2).K_m(1).K_nn;
                J(N1+1:N1+N2,N1+N2+1:N1+N2+C1)        = patches(2).K_m(1).K_nf*c;
                J(N1+1:N1+N2,N1+N2+C1+1:N1+N2+C1+C2)  = patches(2).K_m(1).K_nc + patches(2).K_m(1).K_nf*c;
                %next the combined parts of the system, beginning with the
                %close to the interface part
                J(N1+N2+1:N1+N2+C1,1:N1)                    = patches(1).K_m(1).K_cn;
                J(N1+N2+1:N1+N2+C1,N1+1:N1+N2)              = patches(2).K_m(1).K_cn;
                J(N1+N2+1:N1+N2+C1,N1+N2+1:N1+N2+C1)        = patches(1).K_m(1).K_cc + (patches(1).K_m(1).K_cf+patches(2).K_m(1).K_cf)*c;
                J(N1+N2+1:N1+N2+C1,N1+N2+C1+1:N1+N2+C1+C2)  = patches(2).K_m(1).K_cc + (patches(1).K_m(1).K_cf+patches(2).K_m(1).K_cf)*c;
                %the interface part
                J(N1+N2+C1+1:N1+N2+C1+C2,1:N1)                    = patches(1).K_m(1).K_fn;
                J(N1+N2+C1+1:N1+N2+C1+C2,N1+1:N1+N2)              = patches(2).K_m(1).K_fn;
                J(N1+N2+C1+1:N1+N2+C1+C2,N1+N2+1:N1+N2+C1)        = patches(1).K_m(1).K_fc + (patches(1).K_m(1).K_ff+patches(2).K_m(1).K_ff)*c;
                J(N1+N2+C1+1:N1+N2+C1+C2,N1+N2+C1+1:N1+N2+C1+C2)  = patches(2).K_m(1).K_fc + (patches(1).K_m(1).K_ff+patches(2).K_m(1).K_ff)*c;
                %determine the displacements increments
                du = -J\rsd;
                %write the dislacements inside the patch structure
                patches(1).u(intersect(patches(1).Int_Dofs(1).KP_n,patches(1).freeDofs)) = patches(1).u(intersect(patches(1).Int_Dofs(1).KP_n,patches(1).freeDofs)) + du(1:N1);
                patches(2).u(intersect(patches(2).Int_Dofs(1).KP_n,patches(2).freeDofs)) = patches(2).u(intersect(patches(2).Int_Dofs(1).KP_n,patches(2).freeDofs)) + du(N1+1:N1+N2);
                patches(1).u(intersect(patches(1).Int_Dofs(1).KP_c,patches(1).freeDofs)) = patches(1).u(intersect(patches(1).Int_Dofs(1).KP_c,patches(1).freeDofs)) + du(N1+N2+1:N1+N2+C1);
                patches(2).u(intersect(patches(2).Int_Dofs(1).KP_c,patches(2).freeDofs)) = patches(2).u(intersect(patches(2).Int_Dofs(1).KP_c,patches(2).freeDofs)) + du(N1+N2+C1+1:N1+N2+C1+C2);
                
                patches(1).u(intersect(patches(1).Int_Dofs(1).KP_f,patches(1).freeDofs)) = c* (patches(1).u(intersect(patches(1).Int_Dofs(1).KP_c,patches(1).freeDofs)) + patches(2).u(intersect(patches(2).Int_Dofs(1).KP_c,patches(2).freeDofs)) );
                patches(2).u(intersect(patches(2).Int_Dofs(1).KP_f,patches(2).freeDofs)) = patches(1).u(intersect(patches(1).Int_Dofs(1).KP_f,patches(1).freeDofs));
                
                iter = iter +1;
                if iter == maxiter
                    close all
                    error('No convergence in Newton scheme!')
                end
           end
    else %1 material examples -> C^1 continuity is ensured with the VUKIMS method
        %begin residuum
        %initialise the residuum for all not combined parts of the body
        rsd = [];
        for pl = 1:np
            rsd_p = patches(pl).fint(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs)) - patches(pl).fext(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs));
            rsd = [rsd; rsd_p];
        end
        A = size(rsd,1);
        %now all combined parts of the residuum for every interface
        for pl = 1:size(IM,1)
            ma=IM(pl,4);%master patch
            for i = 1:size(patches(ma).Int,1)
                if pl == patches(ma).Int(i,1)
                    s =  patches(ma).Int(i,5); %slave patch number
                    %position inside the slave patch interface matrix
                    for l = 1:size(patches(s).Int,1)
                        if pl == patches(s).Int(l,1)
                            p = l; %position inside the slave interface matrix
                        else
                        end
                    end
%                     %calculate the extension operators for the current Interface
%                     if D1_flag == 1
%                         Csm = patches(s).Extract_op(p).Csm_1D;
%                         Ts1m1 = pinv(patches(s).Extract_op(p).A1_1D) *patches(ma).Extract_op(i).A1_1D;
%                         Ts2m2 = -pinv(patches(s).Extract_op(p).A2_1D) *Csm*patches(ma).Extract_op(i).A2_1D;
%                     else
%                         Csm = patches(s).Extract_op(p).Csm;
%                         Ts1m1 = pinv(patches(s).Extract_op(p).A1) *patches(ma).Extract_op(i).A1;
%                         Ts2m2 = -pinv(patches(s).Extract_op(p).A2) *Csm*patches(ma).Extract_op(i).A2;
%                     end
                    %first the residuum for the c part
                    rsd_p = patches(ma).fint(intersect(patches(ma).Int_Dofs(i).KP_c,patches(ma).freeDofs)) - patches(ma).fext(intersect(patches(ma).Int_Dofs(i).KP_c,patches(ma).freeDofs)) +...
                        patches(s).fint(intersect(patches(s).Int_Dofs(p).KP_c,patches(s).freeDofs)) - patches(s).fext(intersect(patches(s).Int_Dofs(p).KP_c,patches(s).freeDofs));
                    rsd = [rsd; rsd_p];
                    %second the residuum for the f part
                    rsd_p = patches(s).fint(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs)) - patches(s).fext(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs));
                    
                    rsd_p = patches(ma).fint(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) - patches(ma).fext(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) +...
                        patches(s).fint(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs)) - patches(s).fext(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs));
                    rsd = [rsd; rsd_p];
                end
            end
            B = size(rsd,1);
        end
        %end residuum
        rsn = norm(rsd);
        fprintf(1, ' %4d, residuum_norm= %e\n', iter, rsn);
        if rsn > tolerance
            %begin stiffness matrix
            J1 = zeros(A,size(rsd,1));
            J2 = zeros(B-A,size(rsd,1));
            %first fill the uncombined parts of every patch in ascending order
            A = 0; %counting variable
            %first the K_nn matrices of every patch
            for pl = 1:np
                NP = length(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs));%help variable for the sake of better reading
                J1(1+A:NP+A,1+A:NP+A) = patches(pl).K_m(1).K_nn;
                
                %second the K matrices from all n CP to all other interfaces
                %(K_nci and K_nfi)
                for i = 1:size(patches(pl).Int,1)
                    %first the K_nci matrices
                    %check if the current patch is the master patch of the
                    %interface and get the coloums of this matrix, also get the
                    %position inside the masters interface matrix
                    [alpha,omega,ma,pci] = get_coloum_c(patches(pl).Int(i,1),IM,CM);
                    %calculate the expansion operator Ts2m2
                    if ma==pl %the current patch is the master patch of this interface
                        J1(1+A:NP+A,1+alpha:omega) = patches(pl).K_m(i).K_nc;
                        %second the K_nfi matrices
                        %get the coloums and fill in the matrix
                        [alpha,omega,~,~] = get_coloum_f(patches(pl).Int(i,1),IM);
                        J1(1+A:NP+A,1+alpha:omega) = patches(pl).K_m(i).K_nf;
                    elseif ma~=pl %the current patch is not the master patch of this interface
%                         if D1_flag == 1
%                             Csm = patches(pl).Extract_op(i).Csm_1D;
%                             Ts1m1 = pinv(patches(pl).Extract_op(i).A1_1D) *patches(ma).Extract_op(pci).A1_1D;
%                             Ts2m1 = pinv(patches(pl).Extract_op(i).A2_1D)*(patches(pl).Extract_op(i).A1_1D*pinv(patches(pl).Extract_op(i).A1_1D)+Csm)*patches(ma).Extract_op(pci).A1_1D;
%                             Ts2m2 = -pinv(patches(pl).Extract_op(i).A2_1D) *Csm*patches(ma).Extract_op(pci).A2_1D;
%                         else
%                             Csm = patches(pl).Extract_op(i).Csm;
%                             Ts1m1 = pinv(patches(pl).Extract_op(i).A1) *patches(ma).Extract_op(pci).A1;
%                             Ts2m1 = pinv(patches(pl).Extract_op(i).A2)*(patches(pl).Extract_op(i).A1*pinv(patches(pl).Extract_op(i).A1)+Csm)*patches(ma).Extract_op(pci).A1;
%                             Ts2m2 = -pinv(patches(pl).Extract_op(i).A2) *Csm*patches(ma).Extract_op(pci).A2;
%                         end
                        J1(1+A:NP+A,1+alpha:omega) = patches(pl).K_m(i).K_nc;
                        %second the K_nfi matrices
                        %get the coloums and fill in the matrix
                        [alpha,omega,~,~] = get_coloum_f(patches(pl).Int(i,1),IM);
                        J1(1+A:NP+A,1+alpha:omega) = patches(pl).K_m(i).K_nc + patches(pl).K_m(i).K_nf;
                    end
                    
                    
                end
                A = A + NP;
            end
            %second all the combined parts going over the interface
            %starting with c part of the interface
            A = 0; %counting variable
            for il = 1:size(IM,1)%loop over all interfaces
                m = IM(il,4); %master patch of the current interface
                for i = 1:size(patches(m).Int,1) %finding the interface inside the patch inteface matrix
                    if il == patches(m).Int(i,1)
                        s = patches(m).Int(i,5); %slave patch
                        NF = length(intersect(patches(m).Int_Dofs(i).KP_c,patches(m).freeDofs));%help variable for the sake of better reading
                        %first fill in the two K_cni matrices
                        %first the master K_cni matrix
                        %find the coloumns of the master patch
                        [alpha,omega] = get_coloum_n(PM,m);
                        J2(1+A:NF+A,1+alpha:omega) = patches(m).K_m(i).K_cn;
                        %second the slave K_cni matrix
                        %find the coloumns of the slave patch
                        [alpha,omega] = get_coloum_n(PM,s);
                        %find the position inside the slaves interface matrix
                        for l = 1:size(patches(s).Int,1)
                            if il == patches(s).Int(l,1)
                                p = l; %position inside the slave interface matrix
                            else
                            end
                        end
                        %calculate the extension operators
%                         if D1_flag == 1
%                             Csm = patches(s).Extract_op(p).Csm_1D;
%                             Ts1m1 = pinv(patches(s).Extract_op(p).A1_1D)*patches(m).Extract_op(i).A1_1D;
%                             Ts2m1 = pinv(patches(s).Extract_op(p).A2_1D)*(patches(s).Extract_op(p).A1_1D*pinv(patches(s).Extract_op(p).A1_1D)+Csm)*patches(m).Extract_op(i).A1_1D;
%                             Ts2m2 = -pinv(patches(s).Extract_op(p).A2_1D) *Csm*patches(m).Extract_op(i).A2_1D;
%                         else
%                             Csm = patches(s).Extract_op(p).Csm;
%                             Ts1m1 = pinv(patches(s).Extract_op(p).A1)*patches(m).Extract_op(i).A1;
%                             Ts2m1 = pinv(patches(s).Extract_op(p).A2)*(patches(s).Extract_op(p).A1*pinv(patches(s).Extract_op(p).A1)+Csm)*patches(m).Extract_op(i).A1;
%                             Ts2m2 = -pinv(patches(s).Extract_op(p).A2) *Csm*patches(m).Extract_op(i).A2;
%                         end
                        %fill in the matrix
                        J2(1+A:NF+A,1+alpha:omega) = patches(s).K_m(p).K_cn;
                        %next the combined K_cc matrices on the current interface
                        %find the indicies first
                        [alpha,omega,~,~] = get_coloum_c(il,IM,CM);
                        %fill in the matrices
                        J2(1+A:NF+A,1+alpha:omega) = patches(m).K_m(i).K_cc +  patches(s).K_m(p).K_cc;
                        %after this the K_cf matrices on the current interface
                        %find the indicies first
                        [alpha,omega,~,~] = get_coloum_f(il,IM);
                        %fill in the matrices
                        J2(1+A:NF+A,1+alpha:omega) = patches(m).K_m(i).K_cf +  (patches(s).K_m(p).K_cc+patches(s).K_m(p).K_cf);
                        %now all stiffnes matrices between the master c points
                        %and other c and f points of other interfaces
                        for l = 1:size(patches(m).Int,1)
                            if il ~= patches(m).Int(l,1)
                                %check if the current patch is the master of this
                                %patch
                                [alpha,omega,ma,pci] = get_coloum_c(patches(m).Int(l,1),IM,CM);
                                if ma ==m %the master patch is also the master patch of the other interface
                                    %fill in the matrices
                                    %first the K_cci matrix
                                    K_cci = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_c,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_c,patches(m).freeDofs));
                                    %then the K_cfi matrix
                                    K_cfi = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_c,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_f,patches(m).freeDofs));
                                    %fill in the matrices
                                    %first the K_cci matrices
                                    %coloums already calculated
                                    J2(1+A:NF+A,1+alpha:omega) =  K_cci;
                                    %second the K_cfi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(m).Int(l,1),IM);
                                    %fill in the matrix
                                    J2(1+A:NF+A,1+alpha:omega) = K_cfi;
                                    
                                elseif ma~=m %the master patch is the slave patch of the other interface
                                    %first calculate the extension operator
%                                     if D1_flag == 1
%                                         Csm = patches(m).Extract_op(l).Csm_1D;
%                                         Ts1m1 = pinv(patches(m).Extract_op(l).A1_1D)*patches(ma).Extract_op(pci).A1_1D;
%                                         Ts2m1 = pinv(patches(m).Extract_op(l).A2_1D)*(patches(m).Extract_op(l).A1_1D*pinv(patches(m).Extract_op(l).A1_1D)+Csm)*patches(ma).Extract_op(pci).A1_1D;
%                                         Ts2m2 = -pinv(patches(m).Extract_op(l).A2_1D) *Csm*patches(ma).Extract_op(pci).A2_1D;
%                                     else
%                                         Csm = patches(m).Extract_op(l).Csm;
%                                         Ts1m1 = pinv(patches(m).Extract_op(l).A1)*patches(ma).Extract_op(pci).A1;
%                                         Ts2m1 = pinv(patches(m).Extract_op(l).A2)*(patches(m).Extract_op(l).A1*pinv(patches(m).Extract_op(l).A1)+Csm)*patches(ma).Extract_op(pci).A1;
%                                         Ts2m2 = -pinv(patches(m).Extract_op(l).A2) *Csm*patches(ma).Extract_op(pci).A2;
%                                     end
                                    %fill in the matrices
                                    %first the K_cci matrices
                                    %coloums already calculated
                                    %get the stiffnes matrix
                                    K_cci = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_c,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_c,patches(m).freeDofs));
                                    J2(1+A:NF+A,1+alpha:omega) =  K_cci;
                                    %second the K_cfi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(m).Int(l,1),IM);
                                    %get the stiffnes matrix
                                    K_cfi = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_c,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_f,patches(m).freeDofs));
                                    %fill in the matrix
                                    J2(1+A:NF+A,1+alpha:omega) = K_cci + K_cfi;
                                    
                                end
                                
                            end
                        end
                        %now all stiffnes matrices between the slave c points
                        %and other c and f points of other interfaces
                        %recalculate the extension operator of this interface
%                         if D1_flag == 1
%                             Csm = patches(s).Extract_op(p).Csm_1D;
%                             Ts2m2I = -pinv(patches(s).Extract_op(p).A2_1D) *Csm*patches(m).Extract_op(i).A2_1D;
%                         else
%                             Csm = patches(s).Extract_op(p).Csm;
%                             Ts2m2I = -pinv(patches(s).Extract_op(p).A2) *Csm*patches(m).Extract_op(i).A2;
%                         end
                        for l = 1:size(patches(s).Int,1)
                            if il ~= patches(s).Int(l,1)
                                %check if the current patch is the master of this
                                %interface
                                [alpha,omega,ma,pci] = get_coloum_c(patches(s).Int(l,1),IM,CM);
                                if ma ==s %the slave patch is the master patch of the other interface
                                    %fill in the matrices
                                    %first the K_cci matrix
                                    K_cci = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_c,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_c,patches(s).freeDofs));
                                    %then the K_cfi matrix
                                    K_cfi = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_c,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_f,patches(s).freeDofs));
                                    %fill in the matrices
                                    %first the K_cci matrices
                                    %coloums already calculated
                                    J2(1+A:NF+A,1+alpha:omega) =  K_cci;
                                    %second the K_cfi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(s).Int(l,1),IM);
                                    %fill in the matrix
                                    J2(1+A:NF+A,1+alpha:omega) = K_cfi;
                                    
                                elseif ma~=m %the slave patch is also the slave patch of the other interface
                                    %first calculate the extension operator
%                                     if D1_flag == 1
%                                         Csm = patches(s).Extract_op(l).Csm_1D;
%                                         Ts1m1 = pinv(patches(s).Extract_op(l).A1_1D)*patches(ma).Extract_op(pci).A1_1D;
%                                         Ts2m1 = pinv(patches(s).Extract_op(l).A2_1D)*(patches(s).Extract_op(l).A1_1D*pinv(patches(s).Extract_op(l).A1_1D)+Csm)*patches(ma).Extract_op(pci).A1_1D;
%                                         Ts2m2 = -pinv(patches(s).Extract_op(l).A2_1D) *Csm*patches(ma).Extract_op(pci).A2_1D;
%                                     else
%                                         Csm = patches(s).Extract_op(l).Csm;
%                                         Ts1m1 = pinv(patches(s).Extract_op(l).A1)*patches(ma).Extract_op(pci).A1;
%                                         Ts2m1 = pinv(patches(s).Extract_op(l).A2)*(patches(s).Extract_op(l).A1*pinv(patches(s).Extract_op(l).A1)+Csm)*patches(ma).Extract_op(pci).A1;
%                                         Ts2m2 = -pinv(patches(s).Extract_op(l).A2) *Csm*patches(ma).Extract_op(pci).A2;
%                                     end
                                    %fill in the matrices
                                    %first the K_cci matrices
                                    %coloums already calculated
                                    %get the stiffnes matrix
                                    K_cci = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_c,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_c,patches(s).freeDofs));
                                    J2(1+A:NF+A,1+alpha:omega) =   K_cci;
                                    %second the K_cfi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(s).Int(l,1),IM);
                                    %get the stiffnes matrix
                                    K_cfi = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_c,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_f,patches(s).freeDofs));
                                    %fill in the matrix
                                    J2(1+A:NF+A,1+alpha:omega) =  (K_cci + K_cfi);
                                    
                                end
                            end
                            
                        end
                        %now the f part of the interface
                        %first fill in the two K_fni matrices
                        %first the master K_fni matrix
                        %find the coloumns of the master patch
                        [alpha,omega] = get_coloum_n(PM,m);
                        J2(1+A+NF:2*NF+A,1+alpha:omega) = patches(m).K_m(i).K_fn;
                        %second the slave K_fni matrix
                        %find the coloumns of the slave patch
                        [alpha,omega] = get_coloum_n(PM,s);
                        %find the position inside the slaves interface matrix
                        for l = 1:size(patches(s).Int,1)
                            if il == patches(s).Int(l,1)
                                p = l; %position inside the slave interface matrix
                            else
                            end
                        end
%                         %calculate the extension operators
%                         if D1_flag == 1
%                             Csm = patches(s).Extract_op(p).Csm_1D;
%                             Ts1m1 = pinv(patches(s).Extract_op(p).A1_1D)*patches(m).Extract_op(i).A1_1D;
%                             Ts2m1 = pinv(patches(s).Extract_op(p).A2_1D)*(patches(s).Extract_op(p).A1_1D*pinv(patches(s).Extract_op(p).A1_1D)+Csm)*patches(m).Extract_op(i).A1_1D;
%                             Ts2m2 = -pinv(patches(s).Extract_op(p).A2_1D) *Csm*patches(m).Extract_op(i).A2_1D;
%                         else
%                             Csm = patches(s).Extract_op(p).Csm;
%                             Ts1m1 = pinv(patches(s).Extract_op(p).A1)*patches(m).Extract_op(i).A1;
%                             Ts2m1 = pinv(patches(s).Extract_op(p).A2)*(patches(s).Extract_op(p).A1*pinv(patches(s).Extract_op(p).A1)+Csm)*patches(m).Extract_op(i).A1;
%                             Ts2m2 = -pinv(patches(s).Extract_op(p).A2) *Csm*patches(m).Extract_op(i).A2;
%                         end
                        %fill in the matrix
                        J2(1+A+NF:2*NF+A,1+alpha:omega) = patches(s).K_m(p).K_fn;
                        %next the combined K_ff matrices on the current interface
                        %find the indicies first
                        [alpha,omega,~,~] = get_coloum_f(il,IM);
                        %fill in the matrices
                        J2(1+A+NF:2*NF+A,1+alpha:omega) = patches(m).K_m(i).K_ff +  (patches(s).K_m(p).K_fc + patches(s).K_m(p).K_ff);
                        %after this the K_fc matrices on the current interface
                        %find the indicies first
                        [alpha,omega,~,~] = get_coloum_c(il,IM,CM);
                        %fill in the matrices
                        J2(1+A+NF:2*NF+A,1+alpha:omega) = patches(m).K_m(i).K_fc +  patches(s).K_m(p).K_fc;
                        %now all stiffnes matrices between the master f points
                        %and other c and f points of other interfaces
                        for l = 1:size(patches(m).Int,1)
                            if il ~= patches(m).Int(l,1)
                                %check if the current patch is the master of this
                                %patch
                                [alpha,omega,ma,pci] = get_coloum_c(patches(m).Int(l,1),IM,CM);
                                if ma ==m %the master patch is also the master patch of the other interface
                                    %fill in the matrices
                                    %first the K_fci matrix
                                    K_fci = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_f,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_c,patches(m).freeDofs));
                                    %then the K_ffi matrix
                                    K_ffi = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_f,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_f,patches(m).freeDofs));
                                    %fill in the matrices
                                    %first the K_fci matrices
                                    %coloums already calculated
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) =  K_fci;
                                    %second the K_ffi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(m).Int(l,1),IM);
                                    %fill in the matrix
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) = K_ffi;
                                    
                                elseif ma~=m %the master patch is the slave patch of the other interface
%                                     %first calculate the extension operator
%                                     if D1_flag == 1
%                                         Csm = patches(m).Extract_op(l).Csm_1D;
%                                         Ts1m1 = pinv(patches(m).Extract_op(l).A1_1D)*patches(ma).Extract_op(pci).A1_1D;
%                                         Ts2m1 = pinv(patches(m).Extract_op(l).A2_1D)*(patches(m).Extract_op(l).A1_1D*pinv(patches(m).Extract_op(l).A1_1D)+Csm)*patches(ma).Extract_op(pci).A1_1D;
%                                         Ts2m2 = -pinv(patches(m).Extract_op(l).A2_1D) *Csm*patches(ma).Extract_op(pci).A2_1D;
%                                     else
%                                         Csm = patches(m).Extract_op(l).Csm;
%                                         Ts1m1 = pinv(patches(m).Extract_op(l).A1)*patches(ma).Extract_op(pci).A1;
%                                         Ts2m1 = pinv(patches(m).Extract_op(l).A2)*(patches(m).Extract_op(l).A1*pinv(patches(m).Extract_op(l).A1)+Csm)*patches(ma).Extract_op(pci).A1;
%                                         Ts2m2 = -pinv(patches(m).Extract_op(l).A2) *Csm*patches(ma).Extract_op(pci).A2;
%                                     end
                                    %fill in the matrices
                                    %first the K_cci matrices
                                    %coloums already calculated
                                    %get the stiffnes matrix
                                    K_fci = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_f,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_c,patches(m).freeDofs));
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) =  K_fci;
                                    %second the K_ffi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(m).Int(l,1),IM);
                                    %get the stiffnes matrix
                                    K_ffi = patches(m).K(intersect(patches(m).Int_Dofs(i).KP_f,patches(m).freeDofs),intersect(patches(m).Int_Dofs(l).KP_f,patches(m).freeDofs));
                                    %fill in the matrix
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) = K_fci + K_ffi;
                                    
                                end
                                
                            end
                        end
                        %now all stiffnes matrices between the slave f points
                        %and other c and f points of other interfaces
%                         %recalculate the extension operator of this interface
%                         if D1_flag == 1
%                             Csm = patches(s).Extract_op(p).Csm_1D;
%                             Ts2m2I = -pinv(patches(s).Extract_op(p).A2_1D) *Csm*patches(m).Extract_op(i).A2_1D;
%                         else
%                             Csm = patches(s).Extract_op(p).Csm;
%                             Ts2m2I = -pinv(patches(s).Extract_op(p).A2) *Csm*patches(m).Extract_op(i).A2;
%                         end
                        for l = 1:size(patches(s).Int,1)
                            if il ~= patches(s).Int(l,1)
                                %check if the current patch is the master of this
                                %patch
                                [alpha,omega,ma,pci] = get_coloum_c(patches(s).Int(l,1),IM,CM);
                                if ma ==s %the slave patch is the master patch of the other interface
                                    %fill in the matrices
                                    %first the K_fci matrix
                                    K_fci = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_c,patches(s).freeDofs));
                                    %then the K_ffi matrix
                                    K_ffi = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_f,patches(s).freeDofs));
                                    %fill in the matrices
                                    %first the K_fci matrices
                                    %coloums already calculated
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) =  K_fci;
                                    %second the K_ffi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(s).Int(l,1),IM);
                                    %fill in the matrix
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) = K_ffi;
                                    
                                elseif ma~=s %the slave patch is also the slave patch of the other interface
%                                     %first calculate the extension operators
%                                     if D1_flag == 1
%                                         Csm = patches(s).Extract_op(l).Csm_1D;
%                                         Ts1m1 = pinv(patches(s).Extract_op(l).A1_1D)*patches(ma).Extract_op(pci).A1_1D;
%                                         Ts2m1 = pinv(patches(s).Extract_op(l).A2_1D)*(patches(s).Extract_op(l).A1_1D*pinv(patches(s).Extract_op(l).A1_1D)+Csm)*patches(ma).Extract_op(pci).A1_1D;
%                                         Ts2m2 = -pinv(patches(s).Extract_op(l).A2_1D) *Csm*patches(ma).Extract_op(pci).A2_1D;
%                                     else
%                                         Csm = patches(s).Extract_op(l).Csm;
%                                         Ts1m1 = pinv(patches(s).Extract_op(l).A1)*patches(ma).Extract_op(pci).A1;
%                                         Ts2m1 = pinv(patches(s).Extract_op(l).A2)*(patches(s).Extract_op(l).A1*pinv(patches(s).Extract_op(l).A1)+Csm)*patches(ma).Extract_op(pci).A1;
%                                         Ts2m2 = -pinv(patches(s).Extract_op(l).A2) *Csm*patches(ma).Extract_op(pci).A2;
%                                     end
                                    %fill in the matrices
                                    %first the K_fci matrices
                                    %coloums already calculated
                                    %get the stiffnes matrix
                                    K_fci = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_c,patches(s).freeDofs));
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) =   K_fci;
                                    %second the K_ffi matrix
                                    %get the coloums
                                    [alpha,omega,~,~] = get_coloum_f(patches(s).Int(l,1),IM);
                                    %get the stiffnes matrix
                                    K_ffi = patches(s).K(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs),intersect(patches(s).Int_Dofs(l).KP_f,patches(s).freeDofs));
                                    %fill in the matrix
                                    J2(1+A+NF:2*NF+A,1+alpha:omega) =  (K_fci + K_ffi);
                                end
                            end
                        end
                        
                    end
                end
                A = A+ 2*NF;
            end
            J = [J1; J2];
            %calculate the displacement increments
            du = -J\rsd;
            %write the displacments into the patch structure
            %first the u_n displacements for every patch
            for pl = 1:np
                %get the indicies of the patch
                [alpha,omega] = get_coloum_n(PM,pl);
                patches(pl).u(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs)) = patches(pl).u(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs)) + du(alpha+1:omega);
            end
            for pl = 1:size(IM,1)
                ma = IM(pl,4);
                for i = 1:size(patches(ma).Int,1)
                    if pl==patches(ma).Int(i,1)
                        %get the indicies of the master patch displacements
                        [alpha,omega,m,pci] = get_coloum_f(patches(ma).Int(i,1),IM);
                        %write the master displacements inside the patch
                        %structure
                        patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) = ...
                            patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) + du(alpha+1:omega);
                        %get the indicies of the master patch displacements
                        [alpha,omega,m,pci] = get_coloum_c(patches(ma).Int(i,1),IM,CM);
                        %write the master displacements inside the patch
                        %structure
                        patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_c,patches(ma).freeDofs)) = ...
                            patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_c,patches(ma).freeDofs)) + du(alpha+1:omega);
                        %get the slave patch and the indice inside the slave
                        %interface matrix
                        s = patches(ma).Int(i,5);
                        for l = 1:size(patches(s).Int,1)
                            if pl == patches(s).Int(l,1)
                                p =l;
                            end
                        end
                        %calculate the relation operator T and the slave
%                         %displacements inside the interface
%                         if D1_flag == 1
%                             Csm = patches(s).Extract_op(p).Csm_1D;
%                             Ts1m1 = pinv(patches(s).Extract_op(p).A1_1D)*patches(ma).Extract_op(i).A1_1D;
%                             Ts2m1 = pinv(patches(s).Extract_op(p).A2_1D)*(patches(s).Extract_op(p).A1_1D*pinv(patches(s).Extract_op(p).A1_1D)+Csm)*patches(ma).Extract_op(i).A1_1D;
%                             Ts2m2 = -pinv(patches(s).Extract_op(p).A2_1D) *Csm*patches(ma).Extract_op(i).A2_1D;
%                         else
%                             Csm = patches(s).Extract_op(p).Csm;
%                             Ts1m1 = pinv(patches(s).Extract_op(p).A1)*patches(ma).Extract_op(i).A1;
%                             Ts2m1 = pinv(patches(s).Extract_op(p).A2)*(patches(s).Extract_op(p).A1*pinv(patches(s).Extract_op(p).A1)+Csm)*patches(ma).Extract_op(i).A1;
%                             Ts2m2 = -pinv(patches(s).Extract_op(p).A2) *Csm*patches(ma).Extract_op(i).A2;
%                         end
                        patches(s).u(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs)) = ...
                              patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs));
                        patches(s).u(intersect(patches(s).Int_Dofs(p).KP_c,patches(s).freeDofs)) = ...
                              patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) +...
                            patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_c,patches(ma).freeDofs));
                    end
                end
            end
            iter = iter +1;
            if iter == maxiter
                close all
                error('No convergence in Newton scheme!')
            end
        end
    end
end

%multiple patch C0 continuity
if continuity_flag == 0 && local_refinement_flag == 1
    %begin residuum
    %initialise the residuum for all not combined parts of the body
    rsd = [];
    for pl = 1:np
        rsd_p = patches(pl).fint(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs)) - patches(pl).fext(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs));
        rsd = [rsd; rsd_p];
    end
    A = size(rsd,1);
    %now all combined parts of the residuum for every interface
    for pl = 1:size(IM,1)
        ma=IM(pl,4);%master patch
        for i = 1:size(patches(ma).Int,1)
            if pl == patches(ma).Int(i,1)
                s =  patches(ma).Int(i,5); %slave patch number
                %position inside the slave patch interface matrix
                for l = 1:size(patches(s).Int,1)
                    if pl == patches(s).Int(l,1)
                        p = l; %position inside the slave interface matrix
                    else
                    end
                end
                %calculate the extension operator for the current Interface
                if  D1_flag == 1 %calculation of the relation operator
                    Ts1m1 = pinv(patches(s).Extract_op(p).A1_1D) *patches(ma).Extract_op(i).A1_1D;
                else
                    Ts1m1 = pinv(patches(s).Extract_op(p).A1) *patches(ma).Extract_op(i).A1;
                end
                %set the extension operator to a empty matrix
                Ts1m1 = [];
                rsd_p = patches(ma).fint(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) - patches(ma).fext(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) +...
                    Ts1m1'*patches(s).fint(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs)) - Ts1m1'*patches(s).fext(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs));
                rsd = [rsd; rsd_p];
            end
        end
        B = size(rsd,1);
    end
    %end residuum
    rsn = norm(rsd);
    fprintf(1, ' %4d, residuum_norm= %e\n', iter, rsn);
    if rsn > tolerance
        %begin stiffness matrix
        J1 = zeros(A,size(rsd,1));
        J2 = zeros(B-A,size(rsd,1));
        %first fill the uncombined parts of every patch in ascending order
        A = 0; %counting variable
        %first the K_nn matrices of every patch
        for pl = 1:np
            NP = length(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs));%help variable for the sake of better reading
            J1(1+A:NP+A,1+A:NP+A) = patches(pl).K_m(1).K_nn;
            A = A + NP;
        end
        %second the K_nfi matrices of every patch
        B = 0;
        for pl = 1:np
            for i = 1:size(patches(pl).Int,1)
                NP = length(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs));%help variable for the sake of better reading
                IN = patches(pl).Int(i,1); %number of the interface
                [alpha,omega,m,p] = get_coloum_f(IN,IM); %indicies inside the stiffnes matrix
                if pl == patches(pl).Int(i,3) | patches(pl).h_refinement_flag == 0 %calculation of the relation operator
                    Ts1m1 = eye(length(intersect(patches(pl).Int_Dofs(i).KP_f,patches(pl).freeDofs)));
                else
                    if D1_flag ==1
                        Ts1m1 = pinv(patches(pl).Extract_op(i).A1_1D)*patches(m).Extract_op(p).A1_1D;
                    else
                        Ts1m1 = pinv(patches(pl).Extract_op(i).A1)*patches(m).Extract_op(p).A1;
                    end
                end
                J1(1+B:NP+B,alpha+1:omega) = patches(pl).K_m(i).K_nf*Ts1m1;
            end
            B = B + NP;
        end
        %third all combined patches going over the interfaces
        B = 0;
        for pl = 1:size(IM,1)
            ma=IM(pl,4);%master patch
            for i = 1:size(patches(ma).Int,1)
                if pl == patches(ma).Int(i,1)
                    NF = length(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)); %better reading
                    s = patches(ma).Int(i,5); %slave patch
                    %first the two K_fin matrices
                    %master patch
                    [alpha,omega] = get_coloum_n(PM,ma);
                    J2(1+B:NF+B,alpha+1:omega) = patches(ma).K_m(i).K_fn;
                    %slave patch
                    %first find the index inside the slave interface
                    %matrix
                    for l = 1:size(patches(s).Int,1)
                        if pl == patches(s).Int(l,1)
                            p =l;
                        end
                    end
                    %second get the extraction operator
                    if patches(s).h_refinement_flag ~= 0
                        if D1_flag == 1
                            Ts1m1 = pinv(patches(s).Extract_op(p).A1_1D) *patches(ma).Extract_op(i).A1_1D;
                        else
                            Ts1m1 = pinv(patches(s).Extract_op(p).A1) *patches(ma).Extract_op(i).A1;
                        end
                    else
                        if D1_flag == 1
                            Ts1m1 = eye(length(patches(s).Int_Dofs(p).KP_f)/2);
                        else
                            Ts1m1 = pinv(patches(s).Extract_op(p).A1) *patches(ma).Extract_op(i).A1;
                        end
                    end
                    %sort the slave patch stiffnes matrix inside the
                    %global stiffnes matrix
                    [alpha,omega] = get_coloum_n(PM,s);
                    J2(1+B:NF+B,alpha+1:omega) = Ts1m1'*patches(s).K_m(p).K_fn;
                    %next the combined stiffnes matrices
                    [alpha,omega,~,~] = get_coloum_f(patches(ma).Int(i,1),IM);
                    J2(1+B:NF+B, alpha+1:omega) = patches(ma).K_m(i).K_ff + Ts1m1'* (patches(s).K_m(p).K_ff *Ts1m1);
                    %now all stiffnes matrices between the interface and
                    %other interfaces of the master patch
                    for j = 1:size(patches(ma).Int,1)
                        if i ~=j
                            %find the indicies of the corrsponding
                            %interface matrix
                            [alpha,omega,m,pfj] = get_coloum_f(patches(ma).Int(j,1),IM);
                            %check if the patch is the master patch of
                            %this interface
                            if m~=ma
                                if D1_flag == 1
                                    Ts1m1fj = pinv(patches(ma).Extract_op(j).A1_1D)*patches(m).Extract_op(pfj).A1_1D;
                                else
                                    Ts1m1fj = pinv(patches(ma).Extract_op(j).A1)*patches(m).Extract_op(pfj).A1;
                                end
                            elseif m==ma
                                if D1_flag ==1
                                    Ts1m1fj = eye(length(patches(ma).Int_Dofs(j).KP_f)/2);
                                else
                                    Ts1m1fj = eye(length(patches(ma).Int_Dofs(j).KP_f));
                                end
                            end
                            K_fifj = patches(ma).K(patches(ma).Int_Dofs(i).KP_f,patches(ma).Int_Dofs(j).KP_f);
                            
                            J2(1+B:NF+B, alpha+1:omega) = K_fifj * Ts1m1fj;
                        end
                    end
                    %now the same for the slave patch
                    %recalculate the extansion operator
                    Ts1m1I = pinv(patches(s).Extract_op(p).A1) *patches(ma).Extract_op(i).A1;
                    for j = 1:size(patches(s).Int,1)
                        if j ~= p
                            %find the master patch and the indicies of
                            %this interface
                            [alpha,omega,m,pfj] = get_coloum_f(patches(s).Int(j,1),IM);
                            %check if the patch is the master patch of
                            %this interface
                            if m~=s
                                if D1_flag == 1
                                    Ts1m1fj = pinv(patches(s).Extract_op(j).A1_1D)*patches(m).Extract_op(pfj).A1_1D;
                                else
                                    Ts1m1fj = pinv(patches(s).Extract_op(j).A1)*patches(m).Extract_op(pfj).A1;
                                end
                            elseif m==s
                                if D1_flag == 1
                                    Ts1m1fj = eye(length(patches(s).Int_Dofs(j).KP_f)/2);
                                else
                                    Ts1m1fj = eye(length(patches(s).Int_Dofs(j).KP_f));
                                end
                            end
                            K_fifj = patches(s).K(patches(s).Int_Dofs(p).KP_f,patches(s).Int_Dofs(j).KP_f);
                            J2(1+B:NF+B, alpha+1:omega) = Ts1m1I'*K_fifj * Ts1m1fj;
                        end
                        
                    end
                    B = B + NF;
                end
            end
        end
        J = [J1; J2];
        %calculate the displacement increments
        du = -J\rsd;
        %write the displacments into the patch structure
        %first the u_n displacements for every patch
        for pl = 1:np
            %get the indicies of the patch
            [alpha,omega] = get_coloum_n(PM,pl);
            patches(pl).u(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs)) = patches(pl).u(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs)) + du((alpha+1):omega);
            
            
        end
        %then the u_f displacements in every interface using the T operator
        for pl = 1:size(IM,1)
            ma = IM(pl,4);
            for i = 1:size(patches(ma).Int,1)
                if pl==patches(ma).Int(i,1)
                    %get the indicies of the master patch displacements
                    [alpha,omega,m,p] = get_coloum_f(patches(ma).Int(i,1),IM);
                    %write the master displacements inside the patch
                    %structure
                    patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) = ...
                        patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) + du(alpha+1:omega);
                    %get the slave patch and the indice inside the slave
                    %interface matrix
                    s = patches(ma).Int(i,5);
                    for l = 1:size(patches(s).Int,1)
                        if pl == patches(s).Int(l,1)
                            p =l;
                        end
                    end
                    %calculate the relation operator T and the slave
                    %displacements inside the interface
                    if D1_flag == 1
                        Ts1m1 = pinv(patches(s).Extract_op(p).A1_1D)*patches(ma).Extract_op(i).A1_1D;
                    else
                        Ts1m1 = pinv(patches(s).Extract_op(p).A1)*patches(ma).Extract_op(i).A1;
                    end
                    patches(s).u(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs)) = ...
                        Ts1m1 *  patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs));
                end
            end
        end
        iter = iter +1;
        if iter == maxiter
            close all
            error('No convergence in Newton scheme!')
        end
        
    end
    
end
%multiple patch C0 continuity and no local refinement
if continuity_flag == 0 && local_refinement_flag == 0
    %begin residuum
    %initialise the residuum for all not combined parts of the body
    rsd = [];
    for pl = 1:np
        rsd_p = patches(pl).fint(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs)) - patches(pl).fext(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs));
        rsd = [rsd; rsd_p];
    end
    A = size(rsd,1);
    %now all combined parts of the residuum for every interface
    for pl = 1:size(IM,1)
        ma=IM(pl,4);%master patch
        for i = 1:size(patches(ma).Int,1)
            if pl == patches(ma).Int(i,1)
                s =  patches(ma).Int(i,5); %slave patch number
                %position inside the slave patch interface matrix
                for l = 1:size(patches(s).Int,1)
                    if pl == patches(s).Int(l,1)
                        p = l; %position inside the slave interface matrix
                    else
                    end
                end
                %calculate the extension operator for the current Interface
                %if  D1_flag == 1 %calculation of the relation operator
                %    Ts1m1 = pinv(patches(s).Extract_op(p).A1_1D) *patches(ma).Extract_op(i).A1_1D;
                %else
                %    Ts1m1 = pinv(patches(s).Extract_op(p).A1) *patches(ma).Extract_op(i).A1;
                %end
                %set the extension operator to a empty matrix
                %Ts1m1 = [];
                rsd_p = patches(ma).fint(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) - patches(ma).fext(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) +...
                    patches(s).fint(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs)) - patches(s).fext(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs));
                rsd = [rsd; rsd_p];
            end
        end
        B = size(rsd,1);
    end
    %end residuum
    rsn = norm(rsd);
    fprintf(1, ' %4d, residuum_norm= %e\n', iter, rsn);
    if rsn > tolerance
        %begin stiffness matrix
        J1 = zeros(A,size(rsd,1));
        J2 = zeros(B-A,size(rsd,1));
        %first fill the uncombined parts of every patch in ascending order
        A = 0; %counting variable
        %first the K_nn matrices of every patch
        for pl = 1:np
            NP = length(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs));%help variable for the sake of better reading
            J1(1+A:NP+A,1+A:NP+A) = patches(pl).K_m(1).K_nn;
            A = A + NP;
        end
        %second the K_nfi matrices of every patch
        B = 0;
        for pl = 1:np
            for i = 1:size(patches(pl).Int,1)
                NP = length(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs));%help variable for the sake of better reading
                IN = patches(pl).Int(i,1); %number of the interface
                [alpha,omega,m,p] = get_coloum_f(IN,IM); %indicies inside the stiffnes matrix
%                 if pl == patches(pl).Int(i,3) | patches(pl).h_refinement_flag == 0 %calculation of the relation operator
%                     Ts1m1 = eye(length(intersect(patches(pl).Int_Dofs(i).KP_f,patches(pl).freeDofs)));
%                 else
%                     if D1_flag ==1
%                         Ts1m1 = pinv(patches(pl).Extract_op(i).A1_1D)*patches(m).Extract_op(p).A1_1D;
%                     else
%                         Ts1m1 = pinv(patches(pl).Extract_op(i).A1)*patches(m).Extract_op(p).A1;
%                     end
%                 end
                J1(1+B:NP+B,alpha+1:omega) = patches(pl).K_m(i).K_nf;
            end
            B = B + NP;
        end
        %third all combined patches going over the interfaces
        B = 0;
        for pl = 1:size(IM,1)
            ma=IM(pl,4);%master patch
            for i = 1:size(patches(ma).Int,1)
                if pl == patches(ma).Int(i,1)
                    NF = length(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)); %better reading
                    s = patches(ma).Int(i,5); %slave patch
                    %first the two K_fin matrices
                    %master patch
                    [alpha,omega] = get_coloum_n(PM,ma);
                    J2(1+B:NF+B,alpha+1:omega) = patches(ma).K_m(i).K_fn;
                    %slave patch
                    %first find the index inside the slave interface
                    %matrix
                    for l = 1:size(patches(s).Int,1)
                        if pl == patches(s).Int(l,1)
                            p =l;
                        end
                    end
                    %second get the extraction operator
%                     if patches(s).h_refinement_flag ~= 0
%                         if D1_flag == 1
%                             Ts1m1 = pinv(patches(s).Extract_op(p).A1_1D) *patches(ma).Extract_op(i).A1_1D;
%                         else
%                             Ts1m1 = pinv(patches(s).Extract_op(p).A1) *patches(ma).Extract_op(i).A1;
%                         end
%                     else
%                         if D1_flag == 1
%                             Ts1m1 = eye(length(patches(s).Int_Dofs(p).KP_f)/2);
%                         else
%                             Ts1m1 = pinv(patches(s).Extract_op(p).A1) *patches(ma).Extract_op(i).A1;
%                         end
%                     end
                    %sort the slave patch stiffnes matrix inside the
                    %global stiffnes matrix
                    [alpha,omega] = get_coloum_n(PM,s);
                    J2(1+B:NF+B,alpha+1:omega) = patches(s).K_m(p).K_fn;
                    %next the combined stiffnes matrices
                    [alpha,omega,~,~] = get_coloum_f(patches(ma).Int(i,1),IM);
                    J2(1+B:NF+B, alpha+1:omega) = patches(ma).K_m(i).K_ff + (patches(s).K_m(p).K_ff );
                    %now all stiffnes matrices between the interface and
                    %other interfaces of the master patch
                    for j = 1:size(patches(ma).Int,1)
                        if i ~=j
                            %find the indicies of the corrsponding
                            %interface matrix
                            [alpha,omega,m,pfj] = get_coloum_f(patches(ma).Int(j,1),IM);
                            %check if the patch is the master patch of
                            %this interface
%                             if m~=ma
%                                 if D1_flag == 1
%                                     Ts1m1fj = pinv(patches(ma).Extract_op(j).A1_1D)*patches(m).Extract_op(pfj).A1_1D;
%                                 else
%                                     Ts1m1fj = pinv(patches(ma).Extract_op(j).A1)*patches(m).Extract_op(pfj).A1;
%                                 end
%                             elseif m==ma
%                                 if D1_flag ==1
%                                     Ts1m1fj = eye(length(patches(ma).Int_Dofs(j).KP_f)/2);
%                                 else
%                                     Ts1m1fj = eye(length(patches(ma).Int_Dofs(j).KP_f));
%                                 end
%                             end
                            K_fifj = patches(ma).K(patches(ma).Int_Dofs(i).KP_f,patches(ma).Int_Dofs(j).KP_f);
                            
                            J2(1+B:NF+B, alpha+1:omega) = K_fifj;
                        end
                    end
                    %now the same for the slave patch
                    %recalculate the extansion operator
                    %Ts1m1I = pinv(patches(s).Extract_op(p).A1) *patches(ma).Extract_op(i).A1;
                    for j = 1:size(patches(s).Int,1)
                        if j ~= p
                            %find the master patch and the indicies of
                            %this interface
                            [alpha,omega,m,pfj] = get_coloum_f(patches(s).Int(j,1),IM);
                            %check if the patch is the master patch of
                            %this interface
%                             if m~=s
%                                 if D1_flag == 1
%                                     Ts1m1fj = pinv(patches(s).Extract_op(j).A1_1D)*patches(m).Extract_op(pfj).A1_1D;
%                                 else
%                                     Ts1m1fj = pinv(patches(s).Extract_op(j).A1)*patches(m).Extract_op(pfj).A1;
%                                 end
%                             elseif m==s
%                                 if D1_flag == 1
%                                     Ts1m1fj = eye(length(patches(s).Int_Dofs(j).KP_f)/2);
%                                 else
%                                     Ts1m1fj = eye(length(patches(s).Int_Dofs(j).KP_f));
%                                 end
%                             end
                            K_fifj = patches(s).K(patches(s).Int_Dofs(p).KP_f,patches(s).Int_Dofs(j).KP_f);
                            J2(1+B:NF+B, alpha+1:omega) = K_fifj;
                        end
                        
                    end
                    B = B + NF;
                end
            end
        end
        J = [J1; J2];
        %calculate the displacement increments
        du = -J\rsd;
        %write the displacments into the patch structure
        %first the u_n displacements for every patch
        for pl = 1:np
            %get the indicies of the patch
            [alpha,omega] = get_coloum_n(PM,pl);
            patches(pl).u(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs)) = patches(pl).u(intersect(patches(pl).Int_Dofs(1).KP_n,patches(pl).freeDofs)) + du((alpha+1):omega);
            
            
        end
        %then the u_f displacements in every interface using the T operator
        for pl = 1:size(IM,1)
            ma = IM(pl,4);
            for i = 1:size(patches(ma).Int,1)
                if pl==patches(ma).Int(i,1)
                    %get the indicies of the master patch displacements
                    [alpha,omega,m,p] = get_coloum_f(patches(ma).Int(i,1),IM);
                    %write the master displacements inside the patch
                    %structure
                    patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) = ...
                        patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs)) + du(alpha+1:omega);
                    %get the slave patch and the indice inside the slave
                    %interface matrix
                    s = patches(ma).Int(i,5);
                    for l = 1:size(patches(s).Int,1)
                        if pl == patches(s).Int(l,1)
                            p =l;
                        end
                    end
                    %calculate the relation operator T and the slave
                    %displacements inside the interface
%                     if D1_flag == 1
%                         Ts1m1 = pinv(patches(s).Extract_op(p).A1_1D)*patches(ma).Extract_op(i).A1_1D;
%                     else
%                         Ts1m1 = pinv(patches(s).Extract_op(p).A1)*patches(ma).Extract_op(i).A1;
%                     end
                    patches(s).u(intersect(patches(s).Int_Dofs(p).KP_f,patches(s).freeDofs)) = ...
                        patches(ma).u(intersect(patches(ma).Int_Dofs(i).KP_f,patches(ma).freeDofs));
                end
            end
        end
        iter = iter +1;
        if iter == maxiter
            close all
            error('No convergence in Newton scheme!')
        end
        
    end
    
end
end

%subfunction to find the coloums of the respective dofs inside the global
%system
function [A,O,m,p] = get_coloum_f(IN,IM)
R = find(IM(:,1) == IN);
A = IM(R,6);
O = IM(R,7);
m = IM(R,4);
p = IM(R,5);
end
function [A,O] = get_coloum_n(PM,PN)
R = find(PM(:,1) == PN);
A = PM(R,2);
O = PM(R,3);
end
function [A,O,m,p] = get_coloum_c(IN,IM,CM)
R = find(IM(:,1) == IN);
A = CM(R,2);
O = CM(R,3);
m = IM(R,4);
p = IM(R,5);
end
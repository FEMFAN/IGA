function patches = v_refinement_operators(patches,np)
% ---------------------------------------------------------------------
% Subroutine v_refinement_operators.m
% calculation of the virtuall refienement operators as shown in literature Coox et al <<A robust patch coupling method for NURBS-based
% isogeometric analysis of non-conforming multipatch surfaces>>
%
% Author:           Lennard Langerbein
% Date  :           15.02.2022
%
% Input:    patches                             - structure with patch information
%           np                                  - number of patches
%
% Output:   patches                             - strucure with  the computed virtuall refienement operators
%
%----------------------------------------------------------------------
for i = 1:np
    for j = 1:size(patches(i).Int,1)
        %find the interface in which the patch is a master patch
        if i == patches(i).Int(j,3)
            %slave index
            s = patches(i).Int(j,5); %slave Index
            for l = 1:size(patches(s).Int,1)
                if i == patches(s).Int(l,2)
                    p = l; %position inside the slave Index
                else
                end
            end
            
            
            %check at wich basis function the Interface is
            if patches(i).Int(j,4) == 0
                Xim = patches(i).ETA;
                Xis = patches(s).ETA;
            elseif patches(i).Int(j,4) ==1
                Xim = patches(i).XI;
                Xis = patches(s).XI;
            end
            %find the different knots in both knot vectors
            C1 = setdiff(Xis,Xim);
            C2 = setdiff(Xim,Xis);
            %insert the knots into each other knot vector
            Xim_ex = sort([Xim C1]);
            Xis_ex = sort([Xis C2]);
            
            %calculate the matrix Ai for the master and the slave patch
            %master patch
            if length(C1) > 0
                %calculate the matrix B which is equivalent to the extension
                %operator defined in Hughes et al. <<ISOGEOMETRIC ANALYSIS TOWARD INTEGRATION OF CAD AND FEA>>
                Bm = get_T(Xim,Xim_ex,patches(i).p);
                %calculate the matrix Wm1 and Wm2
                wm1 = [];
                wm2 = [];
                for k = 2:2:length(patches(i).Int_Dofs(j).KP_f)
                    wm1(k/2) = patches(i).w8(patches(i).Int_Dofs(j).KP_f(k)/2);
                    wm2(k/2) = patches(i).w8(patches(i).Int_Dofs(j).KP_c(k)/2);
                end
                Wm1 = zeros(length(wm1));
                Wm2 = zeros(length(wm2));
                for k = 1:length(wm1)
                    for l = 1:length(wm1)
                        if k==l
                            Wm1(k,l) = wm1(k);
                            Wm2(k,l) = wm2(k);
                        end
                    end
                end
                %calculate the matrix Wv
                wvm1 = Bm*wm1';
                wvm2 = Bm*wm2';
                Wvm1 = zeros(length(wvm1));
                Wvm2 = zeros(length(wvm2));
                for k = 1:length(wvm1)
                    for l = 1:length(wvm1)
                        if k==l
                            Wvm1(k,l) = wvm1(k);
                            Wvm2(k,l) = wvm2(k);
                        end
                    end
                end
                %finally calculate the extraction operator
                Am1 = Wvm1^-1 * (Bm* Wm1);
                Am2 = Wvm2^-1 * (Bm* Wm2);
                
            else
                Am1 = eye(size(Xim,2)-patches(i).p-1);
                Am2 = eye(size(Xim,2)-patches(i).p-1);
            end
            %expand the extraction operator to ndm (in this case 2)
            Am1l = zeros(size(Am1,1)*2,size(Am1,2)*2);
            Am2l = zeros(size(Am2,1)*2,size(Am2,2)*2);
            for k = 1:2:size(Am1l,1)
                for l = 1:2:size(Am1l,2)
                    Am1l(k,l) =  Am1(round(k/2),round(l/2));
                    Am1l(k+1,l+1) = Am1(round(k/2),round(l/2));
                    Am2l(k,l) =  Am2(round(k/2),round(l/2));
                    Am2l(k+1,l+1) = Am2(round(k/2),round(l/2));
                end
            end
            
            
            %slave patch
            if length(C2) > 0
                %calculate the matrix B
                Bs = get_T(Xis,Xis_ex,patches(s).p);
                %calculate the matrix Ws1 and Ws2
                ws1 = [];
                ws2 = [];
                for k = 2:2:length(patches(s).Int_Dofs(p).KP_f)
                    ws1(k/2) = patches(s).w8(patches(s).Int_Dofs(p).KP_f(k)/2);
                    ws2(k/2) = patches(s).w8(patches(s).Int_Dofs(p).KP_c(k)/2);
                end
                Ws1 = zeros(length(ws1));
                Ws2 = zeros(length(ws2));
                for k = 1:length(ws1)
                    for l = 1:length(ws1)
                        if k==l
                            Ws1(k,l) = ws1(k);
                            Ws2(k,l) = ws2(k);
                        end
                    end
                end
                %calculate the matrix Wv
                wvs1 = Bs*ws1';
                wvs2 = Bs*ws2';
                Wvs1 = zeros(length(wvs1));
                Wvs2 = zeros(length(wvs2));
                for k = 1:length(wvs1)
                    for l = 1:length(wvs1)
                        if k==l
                            Wvs1(k,l) = wvs1(k);
                            Wvs2(k,l) = wvs2(k);
                        end
                    end
                end
                %finally calculate the extraction operator
                As1 = Wvs1^-1 * (Bs* Ws1);
                As2 = Wvs2^-1 * (Bs* Ws2);
                
            else
                As1 = eye(size(Xis,2)-patches(i).p-1);
                As2 = eye(size(Xis,2)-patches(i).p-1);
            end
            %expand the extraction operator to ndm (in this case 2)
            As1l = zeros(size(As1,1)*2,size(As1,2)*2);
            As2l = zeros(size(As2,1)*2,size(As2,2)*2);
            for k = 1:2:size(As1l,1)
                for l = 1:2:size(As1l,2)
                    As1l(k,l) =  As1(round(k/2),round(l/2));
                    As1l(k+1,l+1) = As1(round(k/2),round(l/2));
                    As2l(k,l) =  As2(round(k/2),round(l/2));
                    As2l(k+1,l+1) = As2(round(k/2),round(l/2));
                end
            end
            
            KPm1 = zeros(length(patches(i).Int_Dofs(j).KP_f)/2,2);
            KPm2 = zeros(length(patches(i).Int_Dofs(j).KP_f)/2,2);
            for d=2:2:length(patches(i).Int_Dofs(j).KP_f)
                KPm1(d/2,:) = patches(i).KP(patches(i).Int_Dofs(j).KP_f(d)/2,:);
                KPm2(d/2,:) = patches(i).KP(patches(i).Int_Dofs(j).KP_c(d)/2,:);
            end
            KPm1_v = Am1 * KPm1;
            KPm2_v = Am2 * KPm2;
            
            KPs1 = zeros(length(patches(s).Int_Dofs(p).KP_f)/2,2);
            KPs2 = zeros(length(patches(s).Int_Dofs(p).KP_f)/2,2);
            for d=2:2:length(patches(s).Int_Dofs(p).KP_f)
                KPs1(d/2,:) = patches(s).KP(patches(s).Int_Dofs(p).KP_f(d)/2,:);
                KPs2(d/2,:) = patches(s).KP(patches(s).Int_Dofs(p).KP_c(d)/2,:);
            end
            KPs1_v = As1 * KPs1;
            KPs2_v = As2 * KPs2;
            %calculate the Csm matrix
            %the Csm matrix is the difference of the master interface CP to the CP
            %close to the interface compared to the same distance at the slave
            %patch
            Csm = zeros(size(KPs1_v,1));
            for d = 1:size(KPs1_v,1)
                Csm(d,d) = sqrt(((KPs1_v(i,1)-KPs2_v(i,1))^2)+(KPs1_v(i,2)-KPs2_v(i,2))^2) / sqrt(((KPm1_v(i,1)-KPm2_v(i,1))^2)+(KPm1_v(i,2)-KPm2_v(i,2))^2);
            end
            
            %expand the Csm matrix to the ndm (in this case 2)
            Csml = zeros(size(Csm,1)*2,size(Csm,2)*2);
            for k = 1:2:size(Csml,1)
                Csml(k,k) =  Csm(round(k/2),round(k/2));
                Csml(k+1,k+1) = Csm(round(k/2),round(k/2));
            end
            
            %lastly write the extraction operators inside the patch structure
            %write the extraction operator for two dimensions inside the patch structure
            patches(i).Extract_op(j).A1 = Am1l;
            patches(i).Extract_op(j).A2 = Am2l;
            patches(s).Extract_op(p).A1 = As1l;
            patches(s).Extract_op(p).A2 = As2l;
            patches(s).Extract_op(p).Csm = Csml;
            patches(i).Extract_op(j).Csm = Csml;
            %write the extraction operator for one dimension inside the patch structure
            patches(i).Extract_op(j).A1_1D = Am1;
            patches(i).Extract_op(j).A2_1D = Am2;
            patches(s).Extract_op(p).A1_1D = As1;
            patches(s).Extract_op(p).A2_1D = As2;
            patches(s).Extract_op(p).Csm_1D = Csm;
            patches(i).Extract_op(j).Csm_1D = Csm;
            
            %check if any interface Dofs are drcl Dofs (the respective Dofs
            %have to be extracted from the the extraction operators)
            %find the respective Dofs
            [val_m,pos_m] = intersect(patches(i).Int_Dofs(j).KP_f,patches(i).drclDofs);
            [val_s,pos_s] = intersect(patches(s).Int_Dofs(p).KP_f,patches(s).drclDofs);
            %modify the extraction operators by deleting the respective
            %row and coloum
            %for the 2-dimensional case
            D1_flag =0;
            if D1_flag ==0
                patches(i).Extract_op(j).A1(pos_s,:) = [];
                patches(i).Extract_op(j).A1(:,pos_m) = [];
                patches(i).Extract_op(j).A2(pos_s,:) = [];
                patches(i).Extract_op(j).A2(:,pos_m) = [];
                patches(s).Extract_op(p).A1(pos_s,:) = [];
                patches(s).Extract_op(p).A1(:,pos_m) = [];
                patches(s).Extract_op(p).A2(pos_s,:) = [];
                patches(s).Extract_op(p).A2(:,pos_m) = [];
                patches(i).Extract_op(j).Csm(pos_s,:) = [];
                patches(i).Extract_op(j).Csm(:,pos_s) = [];
                patches(s).Extract_op(p).Csm(pos_s,:) = [];
                patches(s).Extract_op(p).Csm(:,pos_s) = [];
            elseif D1_flag == 1
                %and for the 1-dimensional case
                patches(i).Extract_op(j).A1_1D(pos_s,:)= [];
                patches(i).Extract_op(j).A1_1D(:,pos_m)= [];
                patches(i).Extract_op(j).A2_1D(pos_s,:)= [];
                patches(i).Extract_op(j).A2_1D(:,pos_m)= [];
                patches(s).Extract_op(p).A1_1D(pos_s,:) = [];
                patches(s).Extract_op(p).A1_1D(:,pos_m) = [];
                patches(s).Extract_op(p).A2_1D(pos_s,:) = [];
                patches(s).Extract_op(p).A2_1D(:,pos_m) = [];
                patches(s).Extract_op(p).Csm_1D(pos_s,:) = [];
                patches(s).Extract_op(p).Csm_1D(:,pos_s) = [];
                patches(i).Extract_op(j).Csm_1D(pos_s,:) = [];
                patches(i).Extract_op(j).Csm_1D(:,pos_s) = [];
            end
        end
    end
end


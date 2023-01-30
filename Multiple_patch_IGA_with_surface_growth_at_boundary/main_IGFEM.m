% ========================================================================
% Isogeometric Finite Elements Analysis
% M. Sc. Carina Witt
% TU Dortmund
% Institute of Mechanics
% ========================================================================
% Isogeometric Finite Element Code for 2d elements
% ========================================================================
% tic
clc;
close all;
clear all;
format short;
% set paths for all subroutines
%set_paths % for Windows
set_paths_2 % for Linux
%% Input data for geometry
% call input data (control points, weights, knot vectors,...)
% Remark: control points have to be sorted by curves in the input file!
%chose input data
%beam_2_materials_input_1;  %beam with a change of material
%beam_growing_top_input_1;   %beam with a growing top
beam_growing_upper_top_input_1;   %beam with a growing upper top
%D1_C0_C1_input_1;          %1-dimensional beam
%plate_with_hole_input_1;   %plate with hole
%truss_input_1;               %truss with no constant cross section
%% Material model
ndm = 2;
%flag indicating two different materials in the patches
% 0: steel in the whole domain
% 1: first patch is steel. Second patch is aluminum
material_flag = 0;
%load the specific material model for every patch
% 1: steel parameters E=2.1e5, v=0.3
% 2: aluminum parameters E=7e4, v = 0.35
for pl = 1:np
    patches(pl).E_4 = Hooke_material(ndm,pl,material_flag);
end
%% Pre-processing

% reshape list of control points and weights for every patch
for pl =1:np
    patches(pl).B = zeros(patches(pl).n,patches(pl).ndm,patches(pl).m);
    patches(pl).w = zeros(patches(pl).n,1,patches(pl).m);
    for i = 0:(patches(pl).m-1)
        patches(pl).B(:,:,i+1) = patches(pl).KP((i*patches(pl).n+1):((i+1)*patches(pl).n),:);
        patches(pl).w(:,1,i+1) = patches(pl).w8((i*patches(pl).n+1):((i+1)*patches(pl).n),1);
    end
end
% choose the desired conituity in the displacment field
% 0 - C^0 continuity
% 1 - C^1 continuity
continuity_flag =1;
% choose if mesh should be modified by knot removal
knot_removal_flag = 0;
remove_knots=[];

% choose if mesh in any patch should be refined by knot insertion
% and define the new knots for every patch as a row vector
% 0 - no refinement
% 3 - knot insertion active for both knot vectors
% warning: problems may occur if multiplicity of knots is changed! (continuity changes)
% -> p-refinement
%index i for patch i
%new_knot = [   new knots Xi patch 1;   new knots Eta patch 1;
%               new knots Xi patch 2;   new knots Eta patch 2;
%at first all patches have to be refined equally
h_refinement_flag = [3;3;3];
new_knots = [   [1:4]*1/5;   [1:4]*1/5;
    [1:4]*1/5;    [1:4]*1/5;
    [1:4]*1/5;    [1:4]*1/5];
%modifing the patch structure with the  user input and saving the knot vectors
%for the case of refinement
for pl=1:pl
    %writing the user input into the patches structure
    patches(pl).h_refinement_flag = h_refinement_flag(pl);
    patches(pl).new_knots_XI = new_knots(pl*2-1,:);
    patches(pl).new_knots_ETA = new_knots(pl*2,:);
end

% construction and plot of control polygon
% reshape of control points and weights
% mesh refinement is obtained in this subroutine if chosen above
% modification of control points, weights and BC description if mesh is refined
for pl=1:np
    [patches(pl).XI,patches(pl).ETA,patches(pl).nkn_XI,patches(pl).nekn_XI,patches(pl).XI_elem,patches(pl).nkn_ETA,patches(pl).nekn_ETA,patches(pl).ETA_elem,patches(pl).KP,patches(pl).w8,patches(pl).B,patches(pl).w,patches(pl).n,patches(pl).m] ...
        = controlpolygon(patches(pl).n,patches(pl).m,patches(pl).p,patches(pl).q,patches(pl).XI,patches(pl).ETA,patches(pl).nkn_XI,patches(pl).nekn_XI,patches(pl).XI_elem,patches(pl).nkn_ETA,patches(pl).nekn_ETA,patches(pl).ETA_elem,patches(pl).ndm,patches(pl).KP,patches(pl).w8,knot_removal_flag,remove_knots,patches(pl).h_refinement_flag,patches(pl).new_knots_XI,patches(pl).new_knots_ETA);
end
%% local refinement
%choose if any patch should be locally refined
% 0 - no refinemnt
% 2 - double all elements in all directions
% 3 - triple all elements in all directions
lr_refinement_flag = [0;0;0];
%save the knot vectors of the refined patch (for the extension operator)
for pl = 1:np
    if lr_refinement_flag(pl)~=0
        patches(pl).Xi_orig = patches(pl).XI;
        patches(pl).Eta_orig = patches(pl).ETA;
    end
    %if locall refinement is used the desired patch is refined again
    if lr_refinement_flag(pl) ~=0
        %insertion of the new knots
        [patches(pl)] = new_knots_lr(patches(pl),lr_refinement_flag(pl));
        %set the refinement flag to 3 and use the refinement function
        patches(pl).h_refinement_flag = 3;
        [patches(pl).XI,patches(pl).ETA,patches(pl).nkn_XI,patches(pl).nekn_XI,patches(pl).XI_elem,patches(pl).nkn_ETA,patches(pl).nekn_ETA,patches(pl).ETA_elem,patches(pl).KP,patches(pl).w8,patches(pl).B,patches(pl).w,patches(pl).n,patches(pl).m] ...
            = controlpolygon(patches(pl).n,patches(pl).m,patches(pl).p,patches(pl).q,patches(pl).XI,patches(pl).ETA,patches(pl).nkn_XI,patches(pl).nekn_XI,patches(pl).XI_elem,patches(pl).nkn_ETA,patches(pl).nekn_ETA,patches(pl).ETA_elem,patches(pl).ndm,patches(pl).KP,patches(pl).w8,knot_removal_flag,remove_knots,patches(pl).h_refinement_flag,patches(pl).new_knots_XI,patches(pl).new_knots_ETA);
    else
        patches(pl).h_refinement_flag = 0;
    end
end
%% Interface finder
%begin interface finder
for i = 1:np
    %find the interfaces inside the interface matrix for the current patch
    for j = 1:size(patches(i).Int,1)
        %find the interfaces in which the current patch is the master patch
        if patches(i).pn == patches(i).Int(j,3)
            %find the index of the interface inside the slave patch
            for l = 1:size(patches(patches(i).Int(j,2)).Int,1)
                if i == patches(patches(i).Int(j,2)).Int(l,2)
                    s = l;
                    break
                else
                end
            end
            [patches(i).Int_Dofs(j).KP_f,patches(i).Int_Dofs(j).KP_c,patches(i).Int_Dofs(j).KP_n,patches(patches(i).Int(j,2)).Int_Dofs(s).KP_f,patches(patches(i).Int(j,2)).Int_Dofs(s).KP_c,patches(patches(i).Int(j,2)).Int_Dofs(s).KP_n] =...
                FindInterface(patches(i).KP,patches(patches(i).Int(j,2)).KP,patches(i).Int(j,:),patches(i).n,patches(i).m,patches(patches(i).Int(j,2)).n,patches(patches(i).Int(j,2)).m,continuity_flag);
        else
        end
    end
end
%end interface finder

%extract the control points which are in two sorting vectors out of the KP_n points
for pl = 1:np
    for i = 1:size(patches(pl).Int_Dofs)
        for j = 1:size(patches(pl).Int_Dofs)
            if continuity_flag ==1
                nc_common = intersect(patches(pl).Int_Dofs(j).KP_c,patches(pl).Int_Dofs(i).KP_n);
                patches(pl).Int_Dofs(i).KP_n = setxor(patches(pl).Int_Dofs(i).KP_n,nc_common);
                nc_common = intersect(patches(pl).Int_Dofs(j).KP_f,patches(pl).Int_Dofs(i).KP_n);
                patches(pl).Int_Dofs(i).KP_n = setxor(patches(pl).Int_Dofs(i).KP_n,nc_common);
            elseif continuity_flag == 0
                nc_common = intersect(patches(pl).Int_Dofs(j).KP_f,patches(pl).Int_Dofs(i).KP_n);
                patches(pl).Int_Dofs(i).KP_n = setxor(patches(pl).Int_Dofs(i).KP_n,nc_common);
            end
        end
    end
end
%end interface finder
%% Plots
% choose if basis functions and NURBS curves shall be plotted
% plots over the full interval of XI and ETA
plot_flag = 0;
if plot_flag == 1
    for pl=1:np
        figure(101)
        Plot_BasisFunctions(patches(pl).XI,patches(pl).nkn_XI,patches(pl).poly_degree,patches(pl).n)
        figure(102)
        Plot_BasisFunctions(patches(pl).ETA,patches(pl).nkn_ETA,patches(pl).poly_degree,patches(pl).m)
        figure(103)
        Plot_NURBS_curves(patches(pl).XI,patches(pl).nkn_XI,patches(pl).n,patches(pl).m,patches(pl).poly_degree,patches(pl).B,patches(pl).w,ndm)
        grid on
        axis('equal')
    end
end
%plot of the mesh
for pl=1:np
    figure(202)
    hold on
    Plot_surface_mesh(patches(pl).XI,patches(pl).ETA,patches(pl).XI_elem,patches(pl).ETA_elem,patches(pl).nkn_XI,patches(pl).nkn_ETA,patches(pl).poly_degree,patches(pl).B,patches(pl).w,patches(pl).ndm);
    ax = gca;
    ax.XAxis.FontSize = 30;
    ax.YAxis.FontSize = 30;
    xlabel('$x$','fontsize',50,'interpreter','latex')
    ylabel('$y$','fontsize',50,'interpreter','latex')
    grid on
    axis('equal')
    hold off
end
%determine the size variables for every patch
for pl = 1:np
    patches(pl).nel = patches(pl).nekn_XI*patches(pl).nekn_ETA; % number of elements in 2D
    patches(pl).nnp = patches(pl).n*patches(pl).m; % number of global basis functions
    patches(pl).nen = (patches(pl).p+1)*(patches(pl).q+1); % number of local basis functions
    patches(pl).nqp = 4; % choose number of Gauss points for numerical integration
    patches(pl).ndf = patches(pl).ndm; % determine number of degrees of freedom here: only purely mechanical problems considered -> ndf=ndm
    patches(pl).tdm = patches(pl).ndm;  % tensor dimension for postprocessing
    patches(pl).u = zeros(patches(pl).nnp*patches(pl).ndm,1); %initialise field variables
end

%% FE-Analysis

% define time step / load step variables
time = 0; % initialise current time
ttime = 1; %end time
dt = 0.2; %time step
nsteps = floor(ttime/dt) + 1; % number of total steps
% define variables for Newton-Raphson loop
maxiter = 10;               % maximum iterations
tolerance     = 1e-7;             % residuum norm tolerance

%growth 
%first row x coordinates of the growth starting point
%second row y coordinates of the growth starting point
%third row x coordinates of the growth ending point
%fourth row y coordinates of the growth ending point
growth = [0    20;
    1.5    1.5;
    0   20;
    2    2];

% loop over time steps
for step = 1:nsteps
    %% Input data for boundary conditions
    
    % call input data (Dirichlet and Neumann BCs, loadcurves,...)
    % Remark: specific format for boundary conditions has to be used!
    
    % direct prescription of BCs to the control points of every patch
    for pl = 1:np
        [fid,patches(pl).u_pre,patches(pl).f_pre,patches(pl).f_b,ttime,dt,loadcurve,D1_flag] = beam_growing_upper_top_input_2(patches(pl).n,patches(pl).m,patches(pl).p,pl,step);
        %% Build connectivities
        % construct IEN and INC array
        % for connectivity of local and global basis functions and NURBS coordinates
        [patches(pl).IEN,patches(pl).INC] = connectivity_arrays(patches(pl).XI,patches(pl).ETA,patches(pl).p,patches(pl).q,patches(pl).n,patches(pl).m,patches(pl).ndm,patches(pl).nel,patches(pl).nen,patches(pl).nnp);
    end
    
    fprintf(1,'time = %8.4e\n',time)
    %determine the limits of the growth area
    if step>1
        xlow = growth(1,1) + (growth(3,1) - growth(1,1))* dt*(step-1)/ttime;
        xhig = growth(1,2) + (growth(3,2) - growth(1,2))* dt*(step-1)/ttime;
        ylow = growth(2,1) + (growth(4,1) - growth(2,1))* dt*(step-1)/ttime;
        yhig = growth(4,2);% + (growth(4,2) - growth(2,2))* dt*step/ttime;
    else
        xlow = growth(1,1);
        xhig = growth(1,2); 
        ylow = growth(2,1); 
        yhig = growth(4,2);
        
    end
    
    % determine Dirichlet and Neumann dofs for the current time step
    % (literature (Cottrell): ID array)
    for pl = 1:np
        [~,patches(pl).drclDofs,patches(pl).freeDofs,~,patches(pl).fpre,patches(pl).u_d] = ID_array(patches(pl).ndm,patches(pl).nnp,patches(pl).INC,patches(pl).u_pre,patches(pl).f_pre,loadcurve,time);
    end
    %calculate the virtuall refinement operator
    [patches] = v_refinement_operators(patches,np);
    %calculate the coloums of the interfaces
    IM = get_coloum_number_f(IM,patches,np,continuity_flag);
    %calculate the coloums of the control points close to the interfaces
    CM = get_coloum_number_c(CM,IM);
    %calculate the coloums of the points inside each patch
    PM = get_coloum_number_n(PM,patches);
    % initialize variables for Newton loop
    rsn  = 1.0;  % residuum norm
    iter = 0;    % current iteration
    
    %apply incrment of Dirichlet BC's in the first itertation
    for pl = 1:np
        patches(pl).u(patches(pl).drclDofs) = patches(pl).u_d;
    end
    % Newton loop
    while rsn > tolerance
        
        % initialize global residuum
        %rsd  = zeros(ndf*nnp,1);
        
        for pl = 1:np
            
            patches(pl).K = zeros(patches(pl).nnp*patches(pl).ndm,patches(pl).nnp*patches(pl).ndm);
            patches(pl).fint = zeros(patches(pl).nnp*patches(pl).ndm,1);
            patches(pl).fext = zeros(patches(pl).nnp*patches(pl).ndm,1);
            patches(pl).fvol = zeros(patches(pl).nnp*patches(pl).ndm,1);
            patches(pl).frea = zeros(patches(pl).nnp*patches(pl).ndm,1);
            
            % initialize arrays for post-processing
            patches(pl).cellstruct = struct('epsilon',zeros(patches(pl).tdm*patches(pl).tdm,patches(pl).nqp),'sigma',zeros(patches(pl).tdm*patches(pl).tdm,patches(pl).nqp));
            patches(pl).celldata   = repmat(patches(pl).cellstruct,patches(pl).nel,1);
            % element loop
            for e=1:patches(pl).nel
                
                %fprintf('Element e=%i\n',e);
                
                % global dofs belonging to element e
                edof = zeros(patches(pl).ndf,patches(pl).nen);
                for ien=1:patches(pl).nen
                    for idf=1:patches(pl).ndf
                        edof(idf,ien) = patches(pl).ndf*patches(pl).IEN(ien,e)+idf-patches(pl).ndf;
                    end
                end
                gdof = reshape(edof,patches(pl).ndf*patches(pl).nen,1);
                
                % displacements for element e
                ue = patches(pl).u(edof);
                
                % initialize element stiffness matrix and element forces
                K_e    = zeros(patches(pl).nen*patches(pl).ndf,patches(pl).nen*patches(pl).ndf);
                fint_e = zeros(patches(pl).nen*patches(pl).ndf,1);
                fvol_e = zeros(patches(pl).nen*patches(pl).ndf,1);
                
                % get NURBS coordinates of the lower left-hand corner of element e
                % i.e. take the highest local basis function number
                % this number can be found in IEN(1,e)
                ni_e = patches(pl).INC(patches(pl).IEN(1,e),1); % xi-direction
                nj_e = patches(pl).INC(patches(pl).IEN(1,e),2); % eta-direction
                
                % check if element has zero measure
                if (patches(pl).XI(ni_e + 1)-patches(pl).XI(ni_e)==0) || (patches(pl).ETA(nj_e + 1)-patches(pl).ETA(nj_e)==0)
                    error('Check element number!')
                end
                
                % get coordinates and weights of Gauss points
                [xi_tilde,w8_tilde] = gaussian_quadrature(patches(pl).nqp,patches(pl).ndm);
                
                % calculate coordinates xi and eta in parameter space
                % of the nqp gauss points
                [xi,dxi_dxitilde] = parameter2gauss(patches(pl).ndm,patches(pl).XI,patches(pl).ETA,xi_tilde,ni_e,nj_e,patches(pl).nqp);
                
                % calculate all BSpline basis functions at the nqp gauss points
                % as well as their first, second and third derivatives
                [N_e,patches(pl).n,basis_functions_XI]      = BSpline_BasisFunctions(patches(pl).XI, xi(1,:),patches(pl).n,patches(pl).nkn_XI, patches(pl).p,patches(pl).nqp);
                [M_e,patches(pl).m,basis_functions_ETA]     = BSpline_BasisFunctions(patches(pl).ETA,xi(2,:),patches(pl).m,patches(pl).nkn_ETA,patches(pl).q,patches(pl).nqp);
                [dN_dxi,  d2N_dxi2,  d3N_dxi3]  = BSpline_BasisFunctions_derivatives(basis_functions_XI, patches(pl).XI, patches(pl).n,patches(pl).nqp,patches(pl).p);
                [dM_deta, d2M_deta2, d3M_deta3] = BSpline_BasisFunctions_derivatives(basis_functions_ETA,patches(pl).ETA,patches(pl).m,patches(pl).nqp,patches(pl).q);
                
                % calculate local NURBS basis functions from the B-Splines
                % as well as their first, second and third derivatives
                [R_e,dR_dxi_e,~,~] = ...
                    local_NURBS_basis_higher_derivatives_2(patches(pl).ndm,patches(pl).nen,N_e,dN_dxi,d2N_dxi2,d3N_dxi3,M_e,dM_deta,d2M_deta2,d3M_deta3,patches(pl).p,patches(pl).q,ni_e,nj_e,patches(pl).w,patches(pl).nqp);
                
                % loop over Gauss points
                for qp=1:patches(pl).nqp
                    
                    % compute cartesian coodinates of quadrature points
                    % for postprocessing
                    [x_xi,dx_dxi] = parameter2physical(patches(pl).tdm,patches(pl).ndm,patches(pl).p,patches(pl).q,ni_e,nj_e,R_e(:,qp),dR_dxi_e(:,:,qp),patches(pl).B);
                    [dxi_dx]      = inverse_cramer(patches(pl).ndm,dx_dxi); % alternative: dxi_dx = inv(dx_dxi);
                    
                    %wie bekomme ich die pysikalischen Koordinaten eines
                    %ganzen elements?
                    if sign(ylow - x_xi(2,1))==-1
                        Death = 1e-6;
                    else
                        Death = 1;
                    end
                    
                    % calculate jacobian and its determinant and inverse
                    % mapping between parent element space and physical space
                    [J,detJ,invJ] = jacobian(dx_dxi, dxi_dxitilde);
                    
                    % check value of the jacobian
                    if detJ <= 0
                        error('Error: det(J) <= 0 in element e= %d\n', e);
                    end
                    
                    % save reference volume blonging to quadrature point q
                    dV = detJ * w8_tilde(1,qp);
                    patches(pl).celldata(e).dV(qp,1) = dV;
                    
                    % derivatives of the NURBS basis functions
                    % w.r.t. physical coordinates
                    [dR_dx] = derivative_dRdx(patches(pl).ndm,patches(pl).nen,dR_dxi_e(:,:,qp),dxi_dx);
                    
                    % calculate the derivatives of the displacement
                    gradu = ue*dR_dx;
                    
                    % calculate the strains for element e
                    eps = zeros(patches(pl).tdm,patches(pl).tdm);
                    eps(1:patches(pl).ndm,1:patches(pl).ndm)  = 1/2 * (gradu + gradu');
                    
                    % calculate stress tensor in current quadrature point
                    % isotropic elasticity with Hooke's law
                    sig = zeros(patches(pl).tdm,patches(pl).tdm);
                    for i = 1:patches(pl).tdm
                        for j = 1:patches(pl).tdm
                            for k = 1:patches(pl).tdm
                                for l = 1:patches(pl).tdm
                                    sig(i,j) = sig(i,j) + patches(pl).E_4(i,j,k,l)*eps(k,l); %stresses should not be manipulated (Death>1e-3)*
                                end
                            end
                        end
                    end
                    
                    % save stress and strain values for postprocessing
                    patches(pl).celldata(e).epsilon(:,qp) = [eps(1,1) eps(1,2) eps(2,1) eps(2,2)];
                    patches(pl).celldata(e).sigma(:,qp)   = [sig(1,1) sig(1,2) sig(2,1) sig(2,2)];
                    
                    % loop over nodes aa
                    for aa=1:patches(pl).nen
                        
                        % calculate internal forces
                        fint_A = zeros(patches(pl).ndm,1);
                        
                        for i=1:patches(pl).ndm
                            for j=1:patches(pl).ndm
                                fint_A(i) = fint_A(i) + dR_dx(aa,j) * sig(i,j) * detJ * w8_tilde(1,qp);
                            end
                        end
                        
                        % assembly of internal force vector for element e
                        fint_e( (aa-1)*patches(pl).ndf+1:aa*patches(pl).ndf,1 ) = fint_e( (aa-1)*patches(pl).ndf+1:aa*patches(pl).ndf,1 ) + Death*fint_A;
                        
                        % calculate volume forces
                        fvol_A = patches(pl).f_b * R_e(aa,qp) * detJ * w8_tilde(1,qp);
                        % assembly of volume force vector for element e
                        fvol_e( (aa-1)*patches(pl).ndf+1:aa*patches(pl).ndf,1 ) = fvol_e( (aa-1)*patches(pl).ndf+1:aa*patches(pl).ndf,1 ) + fvol_A';
                        
                        % loop over nodes bb
                        for bb=1:patches(pl).nen
                            
                            % calculate local stiffness matrix
                            K_AB = zeros(patches(pl).ndm,patches(pl).ndm);
                            for kk=1:patches(pl).ndm
                                for ll=1:patches(pl).ndm
                                    for mm=1:patches(pl).ndm
                                        for nn=1:patches(pl).ndm
                                            K_AB(ll,mm) = K_AB(ll,mm) + dR_dx(aa,kk) * patches(pl).E_4(kk,ll,mm,nn) ...
                                                * dR_dx(bb,nn)' * detJ * w8_tilde(1,qp);
                                        end
                                    end
                                end
                            end
                            
                            % assembly of local stiffness matrix for element e
                            K_e( (aa-1)*patches(pl).ndf+1:aa*patches(pl).ndf,(bb-1)*patches(pl).ndf+1:bb*patches(pl).ndf ) = ...
                                K_e((aa-1)*patches(pl).ndf+1:aa*patches(pl).ndf,(bb-1)*patches(pl).ndf+1:bb*patches(pl).ndf ) + K_AB*Death;
                            
                        end % loop over bb
                        
                    end % loop over aa
                    
                end % gauss point loop
                
                % assembly of fint_e and fvol_e into global vector fint
                patches(pl).fint(gdof)   = patches(pl).fint(gdof)   + fint_e;
                patches(pl).fvol(gdof)   = patches(pl).fvol(gdof)   + fvol_e;
                % assembly of K_e into global stiffness matrix
                patches(pl).K(gdof,gdof) = patches(pl).K(gdof,gdof) + K_e;
                
            end %element loop
            
            
            % calculate external force vector
            patches(pl).fext(patches(pl).freeDofs) = patches(pl).fvol(patches(pl).freeDofs) + patches(pl).fpre(patches(pl).freeDofs);
            
        end %patch loop
        
        %System is now descriped composed into np subdomains and hast to be
        %combined together
        
        %combin system
        
        %reerange the system of equations. Literature: Hughes Cotrell multiple
        %patches revisited
        for pl = 1:np
            for i = 1:size(patches(pl).Int,1)
                patches(pl).K_m(i).K_nn = patches(pl).K(intersect(patches(pl).Int_Dofs(i).KP_n,patches(pl).freeDofs),intersect(patches(pl).Int_Dofs(i).KP_n,patches(pl).freeDofs));
                patches(pl).K_m(i).K_nc = patches(pl).K(intersect(patches(pl).Int_Dofs(i).KP_n,patches(pl).freeDofs),intersect(patches(pl).Int_Dofs(i).KP_c,patches(pl).freeDofs));
                patches(pl).K_m(i).K_nf = patches(pl).K(intersect(patches(pl).Int_Dofs(i).KP_n,patches(pl).freeDofs),intersect(patches(pl).Int_Dofs(i).KP_f,patches(pl).freeDofs));
                patches(pl).K_m(i).K_cn = patches(pl).K(intersect(patches(pl).Int_Dofs(i).KP_c,patches(pl).freeDofs),intersect(patches(pl).Int_Dofs(i).KP_n,patches(pl).freeDofs));
                patches(pl).K_m(i).K_cc = patches(pl).K(intersect(patches(pl).Int_Dofs(i).KP_c,patches(pl).freeDofs),intersect(patches(pl).Int_Dofs(i).KP_c,patches(pl).freeDofs));
                patches(pl).K_m(i).K_cf = patches(pl).K(intersect(patches(pl).Int_Dofs(i).KP_c,patches(pl).freeDofs),intersect(patches(pl).Int_Dofs(i).KP_f,patches(pl).freeDofs));
                patches(pl).K_m(i).K_fn = patches(pl).K(intersect(patches(pl).Int_Dofs(i).KP_f,patches(pl).freeDofs),intersect(patches(pl).Int_Dofs(i).KP_n,patches(pl).freeDofs));
                patches(pl).K_m(i).K_fc = patches(pl).K(intersect(patches(pl).Int_Dofs(i).KP_f,patches(pl).freeDofs),intersect(patches(pl).Int_Dofs(i).KP_c,patches(pl).freeDofs));
                patches(pl).K_m(i).K_ff = patches(pl).K(intersect(patches(pl).Int_Dofs(i).KP_f,patches(pl).freeDofs),intersect(patches(pl).Int_Dofs(i).KP_f,patches(pl).freeDofs));
            end
        end
        
        %calculations of the displacements
        [patches,rsn,iter] = calculator(patches,np,iter,tolerance,maxiter,continuity_flag,IM,PM,CM,D1_flag,material_flag,lr_refinement_flag);
        
        
    end % Newton loop
    
    % compute the remaining entries of the external force vector
    for lp = 1:np
        patches(lp).fext(patches(lp).drclDofs) = patches(lp).K(patches(lp).drclDofs,patches(lp).freeDofs)*patches(lp).u(patches(lp).freeDofs) ...
            + patches(lp).K(patches(lp).drclDofs,patches(lp).drclDofs)*patches(lp).u(patches(lp).drclDofs);
        % compute the reaction forces at the Dirichlet boundary
        patches(lp).frea(patches(lp).drclDofs) = patches(lp).fext(patches(lp).drclDofs) - patches(lp).fvol(patches(lp).drclDofs);
        patches(lp).frea_reshaped  = [patches(lp).frea(1:2:(end-1)),patches(lp).frea(2:2:end)];
        % sum(frea)
    end
    %check for continuity of the displacements
    %[V] = continuity(patches(1).u,patches(2).u,patches(1).KP,patches(2).KP,patches);
    
    %% Postprocessing for the current step
    
    % reshape data
    for lp = 1:np
        patches(lp).u_reshaped  = [patches(lp).u(1:2:(end-1)),patches(lp).u(2:2:end)];
        patches(lp).KP_deformed = patches(lp).u_reshaped + patches(lp).KP  ;
        patches(lp).B_deformed  = zeros(patches(lp).n,patches(lp).ndm,patches(lp).m);
        for i = 0:(patches(lp).m-1)
            patches(lp).B_deformed(:,:,i+1) = patches(lp).KP_deformed((i*patches(lp).n+1):((i+1)*patches(lp).n),:);
        end
        
        % determine initial and current configuration and connectivity list
        [patches(lp).x_initial,patches(lp).conn] = postprocessing(patches(lp).XI,patches(lp).ETA,patches(lp).XI_elem,patches(lp).ETA_elem,patches(lp).nkn_XI,patches(lp).nekn_XI,patches(lp).nkn_ETA,patches(lp).nekn_ETA,patches(pl).poly_degree,patches(lp).B,patches(lp).w,patches(lp).ndm,patches(lp).nel,patches(lp).n,patches(lp).p);
        [patches(lp).x_current,~]    = postprocessing(patches(lp).XI,patches(lp).ETA,patches(lp).XI_elem,patches(lp).ETA_elem,patches(lp).nkn_XI,patches(lp).nekn_XI,patches(lp).nkn_ETA,patches(lp).nekn_ETA,patches(pl).poly_degree,patches(lp).B_deformed,patches(lp).w,patches(lp).ndm,patches(lp).nel,patches(lp).n,patches(lp).p);
        patches(lp).nnp_new          = size(patches(lp).x_current,1);
        patches(lp).u_new            = patches(lp).x_current-patches(lp).x_initial;
    end
    %patches = postprocessing_eps(patches,np);
    create_output_vtk_mp(fid, np, step, time, 2, patches(lp).nnp_new, patches(lp).nel, patches(lp).nqp, patches(lp).x_initial, patches(lp).u_new, patches(lp).conn, patches(lp).celldata,patches);
    
    % update time
    time = step * dt;
    
end % loop over nsteps

for lp = 1:np
    % plot original and deformed control polygon
    figure(200)
    hold on
    Plot_controlpolygon(patches(lp).B,patches(lp).w,patches(lp).n,patches(lp).m,patches(lp).ndm)
    axis('equal')
    hold off
    
    figure(201)
    axis equal
    grid on
    hold on
    Plot_controlpolygon(patches(lp).B_deformed,patches(lp).w,patches(lp).n,patches(lp).m,patches(lp).ndm)
    hold off
    
    % plots of original and deformed surface mesh
    figure(202)
    axis equal
    grid on
    hold on
    Plot_surface_mesh(patches(lp).XI,patches(lp).ETA,patches(lp).XI_elem,patches(lp).ETA_elem,patches(lp).nkn_XI,patches(lp).nkn_ETA,patches(lp).poly_degree,patches(lp).B,patches(lp).w,patches(lp).ndm);
    hold off
    figure(203)
    axis equal
    grid on
    hold on
    Plot_surface_mesh(patches(lp).XI,patches(lp).ETA,patches(lp).XI_elem,patches(lp).ETA_elem,patches(lp).nkn_XI,patches(lp).nkn_ETA,patches(lp).poly_degree,patches(lp).B_deformed,patches(lp).w,patches(lp).ndm);
    hold off
end
%toc
fprintf('IG-Analysis finished\n\n');
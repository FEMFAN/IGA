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

beam_input_1;

%% Material model
ndm = 2;
[E4] = Hooke_material(ndm);

%% Pre-processing

% polynomial order
% here: same polynomial degree in both parametric directions
%p = poly_degree;
%q = p;

% extract elements and number of knots from the knot vectors
%[nkn_XI, nekn_XI,  XI_elem] = element_extraction(XI);
%[nkn_ETA,nekn_ETA,ETA_elem] = element_extraction(ETA);

% number of basis functions
%n = nkn_XI-p-1;
%m = nkn_ETA-q-1;

% reshape list of control points and weights for every patch
for pl =1:np
    patches(pl).B = zeros(patches(pl).n,patches(pl).ndm,patches(pl).m);
    patches(pl).w = zeros(patches(pl).n,1,patches(pl).m);
    for i = 0:(patches(pl).m-1)
        patches(pl).B(:,:,i+1) = patches(pl).KP((i*patches(pl).n+1):((i+1)*patches(pl).n),:);
        patches(pl).w(:,1,i+1) = patches(pl).w8((i*patches(pl).n+1):((i+1)*patches(pl).n),1);
    end
end

% choose if mesh should be modified by knot removal
knot_removal_flag = 0;
remove_knots=[];

% choose if mesh in any patch should be refined by knot insertion
% and define the new knots for every patch as a row vector
% 0 - no refinement
% 1 - knot insertion active for knot vector XI
% 2 - knot insertion active for knot vector ETA
% 3 - knot insertion active for both knot vectors
% warning: problems may occur if multiplicity of knots is changed! (continuity changes)
% -> p-refinement
h_refinement_flag_patch_1 = 0;
h_refinement_flag_patch_2 = 0;

new_knots_XI_patch_1 = [1:29]*1/30;
new_knots_ETA_patch_1 = [1:29]*1/30;
new_knots_XI_patch_2 = [1:29]*1/30;
new_knots_ETA_patch_2 = [1:29]*1/30;

%writing the user input into the patches structure
patches(1).h_refinement_flag = h_refinement_flag_patch_1;
patches(2).h_refinement_flag = h_refinement_flag_patch_2;

patches(1).new_knots_XI = new_knots_XI_patch_1;
patches(1).new_knots_ETA = new_knots_ETA_patch_1;
patches(2).new_knots_XI = new_knots_XI_patch_2;
patches(2).new_knots_ETA = new_knots_ETA_patch_2;

% construction and plot of control polygon
% reshape of control points and weights
% mesh refinement is obtained in this subroutine if chosen above
% modification of control points, weights and BC description if mesh is refined
for pl=1:np
    [patches(pl).XI,patches(pl).ETA,patches(pl).nkn_XI,patches(pl).nekn_XI,patches(pl).XI_elem,patches(pl).nkn_ETA,patches(pl).nekn_ETA,patches(pl).ETA_elem,patches(pl).KP,patches(pl).w8,patches(pl).B,patches(pl).w,patches(pl).n,patches(pl).m] ...
        = controlpolygon(patches(pl).n,patches(pl).m,patches(pl).p,patches(pl).q,patches(pl).XI,patches(pl).ETA,patches(pl).nkn_XI,patches(pl).nekn_XI,patches(pl).XI_elem,patches(pl).nkn_ETA,patches(pl).nekn_ETA,patches(pl).ETA_elem,patches(pl).ndm,patches(pl).KP,patches(pl).w8,knot_removal_flag,remove_knots,patches(pl).h_refinement_flag,patches(pl).new_knots_XI,patches(pl).new_knots_ETA);
end
% Plots
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

for pl=1:np
    figure(202)
    hold on
    Plot_surface_mesh(patches(pl).XI,patches(pl).ETA,patches(pl).XI_elem,patches(pl).ETA_elem,patches(pl).nkn_XI,patches(pl).nkn_ETA,patches(pl).poly_degree,patches(pl).B,patches(pl).w,patches(pl).ndm);
    %grid on
    axis('equal')
    hold off
end
%determine the size variables for every patch
for pl = 1:np
    % number of elements in 2D
    patches(pl).nel = patches(pl).nekn_XI*patches(pl).nekn_ETA;
    % number of global basis functions
    patches(pl).nnp = patches(pl).n*patches(pl).m;
    % number of local basis functions
    patches(pl).nen = (patches(pl).p+1)*(patches(pl).q+1);
    
    % choose number of Gauss points for numerical integration
    patches(pl).nqp = 4;
    
    % determine number of degrees of freedom
    % here: only purely mechanical problems considered -> ndf=ndm
    patches(pl).ndf = patches(pl).ndm;
    
    % tensor dimension for postprocessing
    patches(pl).tdm = patches(pl).ndm;
    
    %initialise field variables
    patches(pl).u = zeros(patches(pl).nnp*patches(pl).ndm,1);
    
    patches(pl).K = zeros(patches(pl).nnp*patches(pl).ndm,patches(pl).nnp*patches(pl).ndm);
    patches(pl).fint = zeros(patches(pl).nnp*patches(pl).ndm,1);
    patches(pl).fext = zeros(patches(pl).nnp*patches(pl).ndm,1);
    patches(pl).fvol = zeros(patches(pl).nnp*patches(pl).ndm,1);
    patches(pl).frea = zeros(patches(pl).nnp*patches(pl).ndm,1);
    
end
%finde the interface
%step 1: find the interface

KP_Interface = intersect(patches(1).KP,patches(2).KP,'rows');

%step 2: find all points on the interface of every single patch

%find the common coordinate (x or y direction)

if KP_Interface(1,1) == KP_Interface(2,1)
    V = KP_Interface(1,1);
    Index = 1;
elseif KP_Interface(1,2) == KP_Interface(2,2)
    V = KP_Interface(1,2);
    Index = 2;
end
patches(1).KP_f = zeros(1);
patches(1).KP_n = zeros(1);
patches(2).KP_f = zeros(1);
patches(2).KP_n = zeros(1);
%find every point in the interface in every patch
for lp = 1:np
    for i = 1:size(patches(lp).KP,1)
        if patches(lp).KP(i,Index) == V
            patches(lp).KP_Interface(end+1,1:2) = patches(lp).KP(i,:);
            patches(lp).KP_Interface(end,3) = i;
            
            if i == 1
                patches(lp).KP_f(end) = i;
                patches(lp).KP_f(end+1)= i + 1;
            else
                if size(patches(lp).KP_f,2) == 1
                    patches(lp).KP_f(end) = i*2 -1;
                    patches(lp).KP_f(end+1) = i*2;
                else
                    
                    patches(lp).KP_f(end+1) = i*2 - 1;
                    patches(lp).KP_f(end+1) = i*2;
                end
            end
        else
            if i == 1
                patches(lp).KP_n(end) = i;
                patches(lp).KP_n(end+1)= i + 1;
            else
                if size(patches(lp).KP_n,2) == 1
                    patches(lp).KP_n(end) = i*2 -1;
                    patches(lp).KP_n(end+1) = i*2;
                else
                    patches(lp).KP_n(end+1) = i*2 - 1;
                    patches(lp).KP_n(end+1) = i*2;
                end
            end
        end
    end
end
for pl = 1:np
    patches(pl).u_n = zeros(length(patches(pl).KP_n),1);
    patches(pl).u_f = zeros(length(patches(pl).KP_f),1);
end

%zeros(patch2.nnp*patch2.ndf,1);

%% Input data for boundary conditions

% call input data (Dirichlet and Neumann BCs, loadcurves,...)
% Remark: specific format for boundary conditions has to be used!

% direct prescription of BCs to the control points of every patch
for pl = 1:np
    [fid,patches(pl).u_pre,patches(pl).f_pre,patches(pl).f_b,ttime,dt,loadcurve] = beam_input_2(patches(pl).n,patches(pl).m,patches(pl).p,pl);
    
    %% Build connectivities
    
    % construct IEN and INC array
    % for connectivity of local and global basis functions and NURBS coordinates
    [patches(pl).IEN,patches(pl).INC] = connectivity_arrays(patches(pl).XI,patches(pl).ETA,patches(pl).p,patches(pl).q,patches(pl).n,patches(pl).m,patches(pl).ndm,patches(pl).nel,patches(pl).nen,patches(pl).nnp);
end
%% FE-Analysis

% define time step / load step variables
time = 0;                     % initialise current time
nsteps = floor(ttime/dt) + 1; % number of total steps
% define variables for Newton-Raphson loop
maxiter = 10;               % maximum iterations
tolerance     = 1e-7;             % residuum norm tolerance

% loop over time steps
for step = 1:nsteps
    
    fprintf(1,'time = %8.4e\n',time)
    
    % determine Dirichlet and Neumann dofs for the current time step
    % (literature (Cottrell): ID array)
    for pl = 1:np
        [~,patches(pl).drclDofs,patches(pl).freeDofs,~,patches(pl).fpre,patches(pl).u_d] = ID_array(patches(pl).ndm,patches(pl).nnp,patches(pl).INC,patches(pl).u_pre,patches(pl).f_pre,loadcurve,time);
    end
    % initialize variables for Newton loop
    rsn  = 1.0;  % residuum norm
    iter = 0;    % current iteration
    
    %apply incrment of Dirichlet BC's in the first itertation
    patches(1).u(patches(1).drclDofs) = patches(1).u_d;
    patches(2).u(patches(2).drclDofs) = patches(2).u_d;
    
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
                
                % fprintf('Element e=%i\n',e);
                
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
                [N_e,n,basis_functions_XI]      = BSpline_BasisFunctions(patches(pl).XI, xi(1,:),patches(pl).n,patches(pl).nkn_XI, patches(pl).p,patches(pl).nqp);
                [M_e,m,basis_functions_ETA]     = BSpline_BasisFunctions(patches(pl).ETA,xi(2,:),patches(pl).m,patches(pl).nkn_ETA,patches(pl).q,patches(pl).nqp);
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
                                    sig(i,j) = sig(i,j) + E4(i,j,k,l)*eps(k,l);
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
                        fint_e( (aa-1)*patches(pl).ndf+1:aa*patches(pl).ndf,1 ) = fint_e( (aa-1)*patches(pl).ndf+1:aa*patches(pl).ndf,1 ) + fint_A;
                        
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
                                            K_AB(ll,mm) = K_AB(ll,mm) + dR_dx(aa,kk) * E4(kk,ll,mm,nn) ...
                                                * dR_dx(bb,nn)' * detJ * w8_tilde(1,qp);
                                        end
                                    end
                                end
                            end
                            
                            % assembly of local stiffness matrix for element e
                            K_e( (aa-1)*patches(pl).ndf+1:aa*patches(pl).ndf,(bb-1)*patches(pl).ndf+1:bb*patches(pl).ndf ) = ...
                                K_e((aa-1)*patches(pl).ndf+1:aa*patches(pl).ndf,(bb-1)*patches(pl).ndf+1:bb*patches(pl).ndf ) + K_AB;
                            
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
        for lp = 1:np
            patches(lp).K_nn = patches(lp).K(intersect(patches(lp).KP_n,patches(lp).freeDofs),intersect(patches(lp).KP_n,patches(lp).freeDofs));
            patches(lp).K_nf = patches(lp).K(intersect(patches(lp).KP_n,patches(lp).freeDofs),intersect(patches(lp).KP_f,patches(lp).freeDofs));
            patches(lp).K_fn = patches(lp).K(intersect(patches(lp).KP_f,patches(lp).freeDofs),intersect(patches(lp).KP_n,patches(lp).freeDofs));
            patches(lp).K_ff = patches(lp).K(intersect(patches(lp).KP_f,patches(lp).freeDofs),intersect(patches(lp).KP_f,patches(lp).freeDofs));
        end
        
        % solve system
        
        %calculate residuum and its norm
        rsd_1 = patches(1).fint(intersect(patches(1).KP_n,patches(1).freeDofs)) - patches(1).fext(intersect(patches(1).KP_n,patches(1).freeDofs));
        rsd_2 = patches(1).fint(intersect(patches(1).KP_f,patches(1).freeDofs)) - patches(1).fext(intersect(patches(1).KP_f,patches(1).freeDofs))...
                + patches(2).fint(intersect(patches(2).KP_f,patches(2).freeDofs)) - patches(2).fext(intersect(patches(2).KP_f,patches(2).freeDofs));
        rsd_3 = patches(2).fint(intersect(patches(2).KP_n,patches(2).freeDofs)) - patches(2).fext(intersect(patches(2).KP_n,patches(2).freeDofs));
        rsd = [rsd_1; rsd_2; rsd_3];
        %modify residuum for inhomogenous Dirichlet Values
        rsn = norm(rsd); %+ norm( patches(1).u_d - patches(1).u(patches(1).drclDofs) ) + norm( patches(2).u_d - patches(2).u(patches(2).drclDofs) ) ;
        fprintf(1, ' %4d, residuum_norm= %e\n', iter, rsn);
        
        %solve system
        if rsn > tolerance
            %calculate increment of node displacements
            %if iter == 0
                %rsd_1 = rsd_1 + patches(1).K(intersect(patches(1).KP_n,patches(1).freeDofs),intersect(patches(1).KP_n,patches(1).drclDofs)) * ...
                %(patches(1).u_d - patches(1).u(patches(1).drclDofs));
                %rsd_2 = rsd_2% +patches(1).K(patches(1).KP_f,patches(1).KP_f) * patches(1).u(patches(1).KP_f) + patches(2).K(patches(2).KP_f,patches(2).KP_f) * patches(2).u(patches(2).KP_f)
                %rsd_2 = rsd_2 + patches(1).K(intersect(patches(1).KP_f,patches(1).freeDofs),intersect(patches(1).KP_f,patches(1).drclDofs)) * ...
                %(patches(1).u_d - patches(1).u(patches(1).drclDofs)) + patches(2).K(intersect(patches(2).KP_f,patches(2).freeDofs),intersect(patches(2).KP_f,patches(2).drclDofs)) * ...
                %(patches(2).u_d - patches(2).u(patches(2).drclDofs));
                %rsd_3 =rsd_3 + patches(2).K(intersect(patches(2).KP_n,patches(2).freeDofs),intersect(patches(2).KP_n,patches(2).drclDofs)) * ...
                %(patches(2).u_d - patches(2).u(patches(2).drclDofs));
               % rhs = [rsd_1; rsd_2; rsd_3];
            %else
                rhs = rsd;
            %end
            %du = zeros(size(rsd_1,1) + size(rsd_2,1) + size(rsd_3,1),1);
            
            null1 = zeros(size(patches(1).K_nf,1),size(patches(1).K_fn,2) + size(patches(1).K_ff,2) + size(patches(2).K_fn,2) - size(patches(1).K_nn,2) - size(patches(1).K_nf,2));
            null2 = zeros(size(patches(2).K_nf,1),size(patches(1).K_fn,2) + size(patches(1).K_ff,2) + size(patches(2).K_fn,2) - size(patches(2).K_nf,2) - size(patches(2).K_nn,2));
            
            J = [patches(1).K_nn patches(1).K_nf null1
                patches(1).K_fn patches(1).K_ff+patches(2).K_ff patches(2).K_fn
                null2 patches(2).K_nf patches(2).K_nn];
            
            %Kspec = Spectraldecomposition(K);
            
            %[U,H,it] = qdwh(K,norm(K,2),min(K),'');
            %Kout,eigvals] = qdwheig(H,norm(H),rank(H),1);
            du = -J \ rhs;
            
            %update displacement values
            %a = patches(1).u_n(intersect(patches(1).KP_n,patches(1).freeDofs))
            %b = du(1:length(intersect(patches(1).KP_n,patches(1).freeDofs)))
            
            patches(1).u(intersect(patches(1).KP_n,patches(1).freeDofs)) = patches(1).u(intersect(patches(1).KP_n,patches(1).freeDofs)) + du(1:length(intersect(patches(1).KP_n,patches(1).freeDofs)));
            patches(1).u(intersect(patches(1).KP_f,patches(1).freeDofs)) = patches(1).u(intersect(patches(1).KP_f,patches(1).freeDofs)) + du(length(intersect(patches(1).KP_n,patches(1).freeDofs)) + 1: (length(intersect(patches(1).KP_n,patches(1).freeDofs)) + length(intersect(patches(1).KP_f,patches(1).freeDofs))) );
            patches(2).u(intersect(patches(2).KP_f,patches(2).freeDofs)) = patches(1).u(intersect(patches(1).KP_f,patches(1).freeDofs));
            patches(2).u(intersect(patches(2).KP_n,patches(2).freeDofs)) = patches(2).u(intersect(patches(2).KP_n,patches(2).freeDofs)) + ...
                du(length(intersect(patches(1).KP_n,patches(1).freeDofs)) + length(intersect(patches(1).KP_f,patches(1).freeDofs)) + 1 : end);
            %patches(1).u_n(intersect(patches(1).KP_n,patches(1).freeDofs)) = patches(1).u_n(intersect(patches(1).KP_n,patches(1).freeDofs)) + du(1:length(intersect(patches(1).KP_n,patches(1).freeDofs)));
            %patches(1).u_f(intersect(patches(1).KP_f,patches(1).freeDofs)) = patches(1).u_f(intersect(patches(1).KP_f,patches(1).freeDofs))' + du(length(intersect(patches(1).KP_n,patches(1).freeDofs)) + 1: (length(intersect(patches(1).KP_n,patches(1).freeDofs)) + length(intersect(patches(1).KP_f,patches(1).freeDofs))) );
            %patches(2).u_f = patches(1).u_f;
            %patches(2).u_n(intersect(patches(2).KP_n,patches(2).freeDofs)) = patches(2).u_n(intersect(patches(2).KP_n,patches(2).freeDofs))' + ...
            %du(length(intersect(patches(1).KP_n,patches(1).freeDofs)) + length(intersect(patches(1).KP_f,patches(1).freeDofs)) + 1 : end);
            %write displacement values at the right spot in the patch
            %displacement vector
            %             for lp = 1:np
            %                 for i = 1:length(patches(lp).u_n)
            %                     if ismember(i, patches(lp).drclDofs)
            %                     else
            %                         patches(lp).u(patches(lp).KP_n(i)) = patches(lp).u_n(i);
            %                     end
            %                 end
            %
            %                 for i = 1:length(patches(lp).u_f)
            %                     if ismember(i, patches(lp).drclDofs)
            %                     else
            %                         patches(lp).u(patches(lp).KP_f(i)) = patches(lp).u_f(i);
            %                     end
            %                 end
            %             end
        end
        % raise the iteration-counter and check if maxiter is reached
        iter = iter +1;
        if iter == maxiter
            close all
            error('No convergence in Newton scheme!')
        end
        
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
    
    
    %% Postprocessing for the current step
    
    % reshape data
    for lp = 1:np
        patches(lp).u_reshaped  = [patches(lp).u(1:2:(end-1)),patches(lp).u(2:2:end)];
        patches(lp).KP_deformed = patches(lp).KP + patches(lp).u_reshaped;
        patches(lp).B_deformed  = zeros(patches(lp).n,patches(lp).ndm,patches(lp).m);
        for i = 0:(patches(lp).m-1)
            patches(lp).B_deformed(:,:,i+1) = patches(lp).KP_deformed((i*patches(lp).n+1):((i+1)*patches(lp).n),:);
        end
        
        % determine initial and current configuration and connectivity list
        [patches(lp).x_initial,patches(lp).conn] = postprocessing(patches(lp).XI,patches(lp).ETA,patches(lp).XI_elem,patches(lp).ETA_elem,patches(lp).nkn_XI,patches(lp).nekn_XI,patches(lp).nkn_ETA,patches(lp).nekn_ETA,patches(pl).poly_degree,patches(lp).B,patches(lp).w,patches(lp).ndm,patches(lp).nel,patches(lp).n,patches(lp).p);
        [patches(lp).x_current,~]    = postprocessing(patches(lp).XI,patches(lp).ETA,patches(lp).XI_elem,patches(lp).ETA_elem,patches(lp).nkn_XI,patches(lp).nekn_XI,patches(lp).nkn_ETA,patches(lp).nekn_ETA,patches(pl).poly_degree,patches(lp).B_deformed,patches(lp).w,patches(lp).ndm,patches(lp).nel,patches(lp).n,patches(lp).p);
        patches(lp).nnp_new          = size(patches(lp).x_current,1);
        patches(lp).u_new            = patches(lp).x_current-patches(lp).x_initial;
        
        % create output files
        %create_output_vtk(fid, step, time, ndf, nnp_new, nel, nqp, x_initial, u_new, conn, celldata);
    end
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
% toc
fprintf('FE-Analysis finished\n\n');

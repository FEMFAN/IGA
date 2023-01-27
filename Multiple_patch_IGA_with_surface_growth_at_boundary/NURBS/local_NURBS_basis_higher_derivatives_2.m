function [R_final,dR_dxi_final,d2R_dxi2_final,d3R_dxi3_final] = local_NURBS_basis_higher_derivatives_2(ndm,nen,N,dN_dxi,d2N_dxi2,d3N_dxi3,M,dM_deta,d2M_deta2,d3M_deta3,p,q,ni,nj,w,nqp)
% ---------------------------------------------------------------------
% Subroutine local_NURBS_basis_higher_derivatives_2.m
% calculates NURBS basis and its first, second and third derivative 
% w.r.t. parametric coordinates xi and eta
%
% Author:           Carina Witt
% Date  :           22.10.2019
%                   28.02.2020: correction for 3rd-order derivatives
%
% Input:    N,M                 - B-Spline basis functions
%           dN_dxi, dM_deta     - B-Spline basis function derivatives
%           d2N_dxi2, d2M_deta2 - B-Spline basis function 2nd order derivatives
%           d3N_dxi3, d3M_deta3 - B-Spline basis function 3rd order derivatives
%           w                   - Weights of the control points
%           ndm                 - number of dimensions
%           nen                 - number of local basis functions
%           p,q                 - polynomial degrees
%           ni,nj               - NURBS coordiantes
%
% Output:   R        - NURBS basis function in 2D
%           dR_dxi   - derivative of NURBS basis function w.r.t. xi 
%           d2R_dxi2 - 2nd order derivative of NURBS basis function w.r.t. xi 
%           d3R_dxi3 - 3rd order derivative of NURBS basis function w.r.t. xi 
%
%-----------------------------------------------------------------------

% literature: Cottrell, (2.30-2.31) and Algorithm 2 
%             The NURBS book, equation (4.19)
R_final           = zeros(nen,nqp);
dR_dxi_final      = zeros(nen,ndm,nqp);
d2R_dxi2_final    = zeros(nen,ndm,ndm,nqp);
d3R_dxi3_final    = zeros(nen,ndm,ndm,ndm,nqp);

% extract active basis functions for the element % 07.12.18
N         = N(ni-p:ni,:);
dN_dxi    = dN_dxi(ni-p:ni,:);
d2N_dxi2  = d2N_dxi2(ni-p:ni,:);
d3N_dxi3  = d3N_dxi3(ni-p:ni,:);
M         = M(nj-q:nj,:);
dM_deta   = dM_deta(nj-q:nj,:);
d2M_deta2 = d2M_deta2(nj-q:nj,:);
d3M_deta3 = d3M_deta3(nj-q:nj,:);

for qp = 1:nqp

    % initialization
    R_tilde      = zeros(nen,1);
    dR_dxi_tilde = zeros(nen,ndm);
    d2R_dxi2_tilde = zeros(nen,ndm); d2R_dxideta_tilde = zeros(nen,ndm);
    d3R_dxi3_tilde = zeros(nen,ndm); d3R_dxi2deta_tilde = zeros(nen,ndm);
    R      = zeros(nen,1);
    dR_dxi = zeros(nen,ndm);
    d2R_dxi2 = zeros(nen,ndm); d2R_dxideta = zeros(nen,ndm);
    d3R_dxi3 = zeros(nen,ndm); d3R_dxi2deta = zeros(nen,ndm);
    W=0; dW_dxi=0; dW_deta=0; 
    d2W_dxi2=0; d2W_deta2=0; d2W_dxideta=0; d2W_detadxi=0;
    d3W_dxi3=0; d3W_deta3=0; d3W_dxi2deta=0; d3W_deta2dxi=0;

    % Rational basis functions
    loc_num=0;
    for j=1:(q+1)
        for i=1:(p+1)
            loc_num = loc_num+1;

            R_tilde(loc_num) = N(p+2-i,qp)*M(q+2-j,qp)*w(ni+1-i,1,nj+1-j); 

            W = W + R_tilde(loc_num);

            dR_dxi_tilde(loc_num,1) = w(ni+1-i,1,nj+1-j)*M(q+2-j,qp)*dN_dxi(p+2-i,qp);
            dR_dxi_tilde(loc_num,2) = w(ni+1-i,1,nj+1-j)*dM_deta(q+2-j,qp)*N(p+2-i,qp);
            
            d2R_dxi2_tilde(loc_num,1) = w(ni+1-i,1,nj+1-j)*M(q+2-j,qp)*d2N_dxi2(p+2-i,qp);
            d2R_dxi2_tilde(loc_num,2) = w(ni+1-i,1,nj+1-j)*d2M_deta2(q+2-j,qp)*N(p+2-i,qp);
            d2R_dxideta_tilde(loc_num,1) = w(ni+1-i,1,nj+1-j)*dM_deta(q+2-j,qp)*dN_dxi(p+2-i,qp);
            d2R_dxideta_tilde(loc_num,2) = w(ni+1-i,1,nj+1-j)*dM_deta(q+2-j,qp)*dN_dxi(p+2-i,qp); %=d2R_dxideta(loc_num,1)
            
            d3R_dxi3_tilde(loc_num,1) = w(ni+1-i,1,nj+1-j)*M(q+2-j,qp)*d3N_dxi3(p+2-i,qp);
            d3R_dxi3_tilde(loc_num,2) = w(ni+1-i,1,nj+1-j)*d3M_deta3(q+2-j,qp)*N(p+2-i,qp);
            d3R_dxi2deta_tilde(loc_num,1) = w(ni+1-i,1,nj+1-j)*dM_deta(q+2-j,qp)*d2N_dxi2(p+2-i,qp);
            d3R_dxi2deta_tilde(loc_num,2) = w(ni+1-i,1,nj+1-j)*d2M_deta2(q+2-j,qp)*dN_dxi(p+2-i,qp); 
            %=d3R_dxideta2, saved as d3R_dxi2deta(loc_num,2)
            
            dW_dxi  = dW_dxi  + dR_dxi_tilde(loc_num,1);
            dW_deta = dW_deta + dR_dxi_tilde(loc_num,2);
            
            d2W_dxi2    = d2W_dxi2    + d2R_dxi2_tilde(loc_num,1);
            d2W_deta2   = d2W_deta2   + d2R_dxi2_tilde(loc_num,2);
            d2W_dxideta = d2W_dxideta + d2R_dxideta_tilde(loc_num,1);
            d2W_detadxi = d2W_detadxi + d2R_dxideta_tilde(loc_num,2); %=d2W_dxideta
            
            d3W_dxi3     = d3W_dxi3    + d3R_dxi3_tilde(loc_num,1);
            d3W_deta3    = d3W_deta3   + d3R_dxi3_tilde(loc_num,2);
            d3W_dxi2deta = d3W_dxi2deta + d3R_dxi2deta_tilde(loc_num,1);
            d3W_deta2dxi = d3W_deta2dxi + d3R_dxi2deta_tilde(loc_num,2);
        end
    end

    for loc_num=1:nen        
        % attention: R_tilde=N*M*w_ij and R=R_tilde/W
        
        R(loc_num)        = R_tilde(loc_num)/W;
        
        dR_dxi(loc_num,1) = ( dR_dxi_tilde(loc_num,1) - R(loc_num) * dW_dxi  ) / W;    
        dR_dxi(loc_num,2) = ( dR_dxi_tilde(loc_num,2) - R(loc_num) * dW_deta ) / W;  
        
        d2R_dxi2(loc_num,1)    = ( d2R_dxi2_tilde(loc_num,1) - 2*dR_dxi(loc_num,1)*dW_dxi   - R(loc_num) * d2W_dxi2 )  / W;  
        d2R_dxi2(loc_num,2)    = ( d2R_dxi2_tilde(loc_num,2) - 2*dR_dxi(loc_num,2)*dW_deta  - R(loc_num) * d2W_deta2 ) / W;           
        d2R_dxideta(loc_num,1) = ( d2R_dxideta_tilde(loc_num,1) - dR_dxi(loc_num,1)*dW_deta -  dR_dxi(loc_num,2)*dW_dxi - R(loc_num)*d2W_dxideta ) / W;     
        d2R_dxideta(loc_num,2) = ( d2R_dxideta_tilde(loc_num,2) - dR_dxi(loc_num,1)*dW_deta -  dR_dxi(loc_num,2)*dW_dxi - R(loc_num)*d2W_detadxi ) / W;  %=d2R_dxideta(loc_num,1)

        d3R_dxi3(loc_num,1)     = ( d3R_dxi3_tilde(loc_num,1) - 3*d2R_dxi2(loc_num,1)*dW_dxi  - 3*dR_dxi(loc_num,1)*d2W_dxi2  - R(loc_num)*d3W_dxi3  ) / W; %=d3R_dxi3
        d3R_dxi3(loc_num,2)     = ( d3R_dxi3_tilde(loc_num,2) - 3*d2R_dxi2(loc_num,2)*dW_deta - 3*dR_dxi(loc_num,2)*d2W_deta2 - R(loc_num)*d3W_deta3 ) / W; %=d3R_deta3     
        d3R_dxi2deta(loc_num,1) = ( d3R_dxi2deta_tilde(loc_num,1) - 2*d2R_dxideta(loc_num,1)*dW_dxi  - dR_dxi(loc_num,2)*d2W_dxi2 ...
                                    - d2R_dxi2(loc_num,1)*dW_deta - 2*dR_dxi(loc_num,1)*d2W_dxideta  - R(loc_num)*d3W_dxi2deta ) / W; %=d3R_dxi2deta 
        d3R_dxi2deta(loc_num,2) = ( d3R_dxi2deta_tilde(loc_num,2) - 2*d2R_dxideta(loc_num,2)*dW_deta - dR_dxi(loc_num,1)*d2W_deta2 ...
                                    - d2R_dxi2(loc_num,2)*dW_dxi  - 2*dR_dxi(loc_num,2)*d2W_detadxi  - R(loc_num)*d3W_deta2dxi ) / W; %=d3R_dxideta2
    end
    % R

% assembly

R_final(:,qp)              = R;                   % R
dR_dxi_final(:,:,qp)       = dR_dxi;              % dR_dxi, dR_deta

d2R_dxi2_final(:,1,1,qp)   = d2R_dxi2(:,1);       % d2R_dxi2
d2R_dxi2_final(:,1,2,qp)   = d2R_dxideta(:,1);    % d2R_dxideta
d2R_dxi2_final(:,2,1,qp)   = d2R_dxideta(:,2);    % d2R_deta_dxi
d2R_dxi2_final(:,2,2,qp)   = d2R_dxi2(:,2);       % d2R_deta2

d3R_dxi3_final(:,1,1,1,qp) = d3R_dxi3(:,1);        % d3R_dxi3
d3R_dxi3_final(:,1,1,2,qp) = d3R_dxi2deta(:,1);    % d3R_dxi2deta 
d3R_dxi3_final(:,1,2,1,qp) = d3R_dxi2deta(:,1);    % d3R_dxidetadxi 
d3R_dxi3_final(:,2,1,1,qp) = d3R_dxi2deta(:,1);    % d3R_detadxi2 
d3R_dxi3_final(:,1,2,2,qp) = d3R_dxi2deta(:,2);    % d3R_dxideta2 
d3R_dxi3_final(:,2,1,2,qp) = d3R_dxi2deta(:,2);    % d3R_detadxideta 
d3R_dxi3_final(:,2,2,1,qp) = d3R_dxi2deta(:,2);    % d3R_deta2dxi
d3R_dxi3_final(:,2,2,2,qp) = d3R_dxi3(:,2);        % d3R_deta3

end %loop over nqp

end % function

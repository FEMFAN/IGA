function [xi,dxi_dxitilde] = parameter2gauss(ndm,XI,ETA,xi_tilde,ni,nj,nqp)
% ---------------------------------------------------------------------
% Subroutine connectivity_parameter2gauss.m
% mapping from parent element to parameter space
%
% Author:           Carina Witt
% Date  :           13.07.2018
%
% Input:    ndm                          - number of dimensions
%           XI,ETA                       - knot vectors
%           xi_tilde,eta_tilde           - coordinates in parent element
%           ni,nj                        - NURBS coordinates
%           nqp                          - number of quadrature points
%
% Output:   xi                           - coordinate vector in parameter space 
%           dxi_dxitilde                 - gradient of mapping from 
%                                           parameter space to parent element
%
%---------------------------------------------------------------------- 

% mapping from parent element to parameter space
% literature: Cottrell (3.A.3 - 3.A.4)
xi      = zeros(ndm,nqp);
xi(1,:) = ((XI(ni+1)-XI(ni))*xi_tilde(1,:)...
         + (XI(ni+1)+XI(ni))) / 2;
xi(2,:) = ((ETA(nj+1)-ETA(nj))*xi_tilde(2,:)...
         + (ETA(nj+1)+ETA(nj))) / 2;
  
% gradient of the mapping 
% literature: Cottrell (Algorithm 3)
dxi_dxitilde      = zeros(ndm,ndm);
% only the diagonal elements are not zero
dxi_dxitilde(1,1) = (XI(ni+1)-XI(ni))/2;
dxi_dxitilde(2,2) = (ETA(nj+1)-ETA(nj))/2;

end %function

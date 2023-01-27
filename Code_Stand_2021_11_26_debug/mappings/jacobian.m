function [J,detJ,invJ] = jacobian(dx_dxi,dxi_dxitilde)
% ---------------------------------------------------------------------
% Subroutine jacobian.m
% calculates mapping dx/dxi_tilde i.e. mapping from physical space
% to parent element
%
% Author:           Carina Witt
% Date  :           18.07.2018
%
% Input:    dx_dxi         - mapping from physical to parameter space
%           dxi_dxitilde   - gradient of mapping from parameter space to parent element
%
% Output:   J              - mapping from physical space to parent element
%           detJ, invJ     - determinant and inverse of jacobian
%
%----------------------------------------------------------------------

% calculation of jacobian
J = dx_dxi * dxi_dxitilde;  

% determinant of Jacobian
detJ = det(J);

% inverse of Jacobian
invJ = inv(J);

end % function
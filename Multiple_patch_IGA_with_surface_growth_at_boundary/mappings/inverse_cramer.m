function [dxi_dx] = inverse_cramer(ndm,dx_dxi)
% ---------------------------------------------------------------------
% Subroutine inverse_cramer
% calculates the inverse mapping dxi_dx out of dx_dxi
%
% Author:           Carina Witt
% Date  :           26.07.2018
%
% Input:    dx_dxi          - mapping from physical to parameter space
%
% Output:   dxi_dx          - inverse mapping

%----------------------------------------------------------------------

% special cramer rule for 2x2 matrix
if ndm == 2
    dxi_dx = 1/( dx_dxi(1,1)*dx_dxi(2,2)-dx_dxi(1,2)*dx_dxi(2,1) ) ...
             * [dx_dxi(2,2) , -dx_dxi(1,2) ; 
                -dx_dxi(2,1) , dx_dxi(1,1)]; 
else 
    error('3 dimensional problem not yet implemented.')
end
% dxi_dx

end % function
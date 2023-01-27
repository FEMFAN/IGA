function [dR_dx] = derivative_dRdx(ndm,nen,dR_dxi,dxi_dx)
% ---------------------------------------------------------------------
% Subroutine derivative_dRdx
% derivative of NURBS basis w.r.t. physical coordinate x
%
% Author:           Carina Witt
% Date  :           18.07.2018
%
% Input:    ndm             - number of dimensions
%           nen             - number of local basis functions
%           dR_dxi          - gradient of NURBS basis functions w.r.t. xi
%           dxi_dx          - mapping from parameter to pysical space
%
% Output:   dR_dx           - gradient of NURBS basis functions w.r.t. x

%----------------------------------------------------------------------
    
% initialization
dR_dx = zeros(nen,ndm);

% derivation
% literature: Cottrell, Algorithm 3
for loc_num = 1:nen
    for i=1:ndm
        for j=1:ndm
            dR_dx(loc_num,i) = dR_dx(loc_num,i) + dR_dxi(loc_num,j)*dxi_dx(j,i);
        end
    end
end
% dR_dx

end % function
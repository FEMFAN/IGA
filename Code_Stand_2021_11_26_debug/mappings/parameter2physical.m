function [x_xi,dx_dxi] = parameter2physical(tdm,ndm,p,q,ni,nj,R,dR_dxi,B)
% ---------------------------------------------------------------------
% Subroutine connectivity_parameter2physical
% mapping from physical to parameter space
%
% Author:           Carina Witt
% Date  :           18.07.2018
%
% Input:    ndm             - number of dimensions
%           p,q             - polynomial degrees
%           ni,nj           - NURBS coordinates
%           dR_dxi          - gradient of NURBS basis functions
%           B               - control points
%
% Output:   dx_dxi          - mapping from physical to parameter space
%           x_xi            - coordinates x(xi)

%----------------------------------------------------------------------
    
% literature: see Cottrell (3.2) + classic FEM

% initialization
x_xi   = zeros(tdm,1);
dx_dxi = zeros(ndm,ndm);

% coordinate transformation
loc_num = 0;
for j=1:(q+1)
    for i=1:(p+1)
        loc_num = loc_num+1;
        for k = 1:ndm
            x_xi(k) = x_xi(k) + B(ni+1-i,k,nj+1-j)*R(loc_num);
        end
    end
end

% derivation
loc_num = 0;
for j=1:(q+1)
    for i=1:(p+1)
        loc_num = loc_num+1;
        for k = 1:ndm
            for l = 1:ndm
                dx_dxi(k,l) = dx_dxi(k,l) + B(ni+1-i,k,nj+1-j)*dR_dxi(loc_num,l);
            end
        end
    end
end

%dx_dxi

end % function
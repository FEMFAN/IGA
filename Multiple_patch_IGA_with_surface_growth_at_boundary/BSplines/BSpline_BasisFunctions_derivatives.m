function [dN_dxi_final, d2N_dxi2_final, d3N_dxi3_final] = BSpline_BasisFunctions_derivatives(basis_functions,XI,n,nqp,p)
% ---------------------------------------------------------------------
% Subroutine BSpline_BasisFunctions_derivatives.m
% calculation of basis function derivatives at location xi
% derivatives up to 3rd order
%
% Author:           Carina Witt
% Date  :           23.08.2019, 22.10.2019 (3rd derivatives)
%
% Input:    basis_functions - struct with basis functions
%           XI              - knot vector
%           n               - number of basis functions in xi-direction
%           nqp             - number of quadrature points
%           p               - polynomial degree
%
% Output:   dN_dxi_final    - first derivative of NURBS basis functions
%           d2N_dxi2_final  - second derivative of NURBS basis functions
%           d3N_dxi3_final  - third derivative of NURBS basis functions
%
%----------------------------------------------------------------------

dN_dxi_final   = zeros(n,nqp);
d2N_dxi2_final = zeros(n,nqp);
d3N_dxi3_final = zeros(n,nqp);

for qp = 1:nqp

    %% calculation of the first derivative w.r.t. xi
    %...%

    dN_dxi = zeros(n,1);
    d2N_dxi2 = zeros(n,1);
    d3N_dxi3 = zeros(n,1);

    for i = 1:n

        % calculate constants a
        a = zeros(4,4); % for first, second and third derivatives k={1,2,3}
        
        for k = 1:2 % for first, second and third derivatives

            l = k+1;

            a(1,1) = 1;
            
            diff1  = XI(i+p+1)-XI(i+k);
            if diff1 ~= 0
                a(l,l) = -a(l-1,l-1)/diff1;
            else 
                a(l,l) = 0;
            end %if
            
            diff2  = XI(i+p-k+1)-XI(i);
            if diff2 ~= 0
                a(l,1) = a(l-1,1)/diff2;
            else
                a(l,1) = 0;
            end %if
            
            for j=1:k-1
                diff3 = XI(i+p+j-k+1)-XI(i+j);
                if diff3 ~= 0
                    a(l,j+1) = (a(l-1,j+1)-a(l-1,j))/diff3;
                else 
                    a(l,j+1) = 0;
                end
            end
        
        end % for
                
        % calculate first derivative
        k=1;
        if k<=p
            for j=0:k
                p_k = p-k;
                if p_k == 0
                    dN_dxi(i,1) = dN_dxi(i,1) + a(k+1,j+1) * basis_functions.N_0(i+j,qp);
                elseif p_k == 1
                    dN_dxi(i,1) = dN_dxi(i,1) + a(k+1,j+1) * basis_functions.N_1(i+j,qp);
                elseif p_k == 2
                    dN_dxi(i,1) = dN_dxi(i,1) + a(k+1,j+1) * basis_functions.N_2(i+j,qp);
                elseif p_k == 3
                    dN_dxi(i,1) = dN_dxi(i,1) + a(k+1,j+1) * basis_functions.N_3(i+j,qp);
                elseif p_k == 4
                    dN_dxi(i,1) = dN_dxi(i,1) + a(k+1,j+1) * basis_functions.N_4(i+j,qp);
                else
                    error('Check polynomial degree')
                end
            end
            dN_dxi(i,1) = factorial(p)/factorial(p-k) * dN_dxi(i,1);
        else 
            dN_dxi(i,1) = 0;
        end % if
        
        % calculate second derivative
        k=2;
        if k<=p
            for j=0:k
                p_k = p-k;
                if p_k == 0
                    d2N_dxi2(i,1) = d2N_dxi2(i,1) + a(k+1,j+1) * basis_functions.N_0(i+j,qp);
                elseif p_k == 1
                    d2N_dxi2(i,1) = d2N_dxi2(i,1) + a(k+1,j+1) * basis_functions.N_1(i+j,qp);
                elseif p_k == 2
                    d2N_dxi2(i,1) = d2N_dxi2(i,1) + a(k+1,j+1) * basis_functions.N_2(i+j,qp);
                elseif p_k == 3
                    d2N_dxi2(i,1) = d2N_dxi2(i,1) + a(k+1,j+1) * basis_functions.N_3(i+j,qp);
                elseif p_k == 4
                    d2N_dxi2(i,1) = d2N_dxi2(i,1) + a(k+1,j+1) * basis_functions.N_4(i+j,qp);                    
                else
                    error('Check polynomial degree')
                end
            end
            d2N_dxi2(i,1) = factorial(p)/factorial(p-k) * d2N_dxi2(i,1);
        else 
            d2N_dxi2(i,1) = 0;
        end % if
        
        % calculate third derivative
        k=3;
        if k<=p
            for j=0:k
                p_k = p-k;
                if p_k == 0
                    d3N_dxi3(i,1) = d3N_dxi3(i,1) + a(k+1,j+1) * basis_functions.N_0(i+j,qp);
                elseif p_k == 1
                    d3N_dxi3(i,1) = d3N_dxi3(i,1) + a(k+1,j+1) * basis_functions.N_1(i+j,qp);
                elseif p_k == 2
                    d3N_dxi3(i,1) = d3N_dxi3(i,1) + a(k+1,j+1) * basis_functions.N_2(i+j,qp);
                elseif p_k == 3
                    d3N_dxi3(i,1) = d3N_dxi3(i,1) + a(k+1,j+1) * basis_functions.N_3(i+j,qp);
                elseif p_k == 4
                    d3N_dxi3(i,1) = d3N_dxi3(i,1) + a(k+1,j+1) * basis_functions.N_4(i+j,qp);                    
                else
                    error('Check polynomial degree')
                end
            end
            d3N_dxi3(i,1) = factorial(p)/factorial(p-k) * d3N_dxi3(i,1);
        else 
            d3N_dxi3(i,1) = 0;
        end % if
    
    end % for

    %assembly
    dN_dxi_final(:,qp)   = dN_dxi;
    d2N_dxi2_final(:,qp) = d2N_dxi2;
    d3N_dxi3_final(:,qp) = d3N_dxi3;
        
end % for

end % function
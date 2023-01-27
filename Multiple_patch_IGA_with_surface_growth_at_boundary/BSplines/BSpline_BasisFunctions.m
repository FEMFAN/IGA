function [N_final,n,basis_functions] = BSpline_BasisFunctions(XI,xi,n,nkn,poly_degree,nqp)
% ---------------------------------------------------------------------
% Subroutine BSpline_BasisFunctions.m
% calculation of basis functions at location xi
% for a polynomial degree up to 4
%
% Author:           Carina Witt
% Date  :           23.08.2019
%
% Input:    xI              - coordinate xi
%           XI              - knot vector
%           n               - number of basis functions in xi-direction
%           nkn             - number of knots in the knot vector XI
%           poly_degree     - polynomial degree of the basis functions
%           nqp             - number of quadrature points
%
% Output:   N               - value of basis functions at location xi
%           n               - number of Basis functions
%           basis_functions - struct with basis functions
%
%----------------------------------------------------------------------
    
%%% recursive calculation of the basis functions  %%%
%%%  and calculation of their derivatives for p<0 %%%
% literature: Cotrell, (2.1 - 2.2) 
% and NURBS book (2.9)

%...
if poly_degree == 4
    basis_functions = struct('N_0',zeros(n+poly_degree,nqp), ...
                             'N_1',zeros(n+poly_degree-1,nqp), ...
                             'N_2',zeros(n+poly_degree-2,nqp), ...
                             'N_3',zeros(n+poly_degree-3,nqp), ...
                             'N_4',zeros(n+poly_degree-4,nqp) );    
elseif poly_degree == 3 
    basis_functions = struct('N_0',zeros(n+poly_degree,nqp), ...
                             'N_1',zeros(n+poly_degree-1,nqp), ...
                             'N_2',zeros(n+poly_degree-2,nqp), ...
                             'N_3',zeros(n+poly_degree-3,nqp) );
elseif poly_degree == 2
        basis_functions = struct('N_0',zeros(n+poly_degree,nqp), ...
                                 'N_1',zeros(n+poly_degree-1,nqp), ...
                                 'N_2',zeros(n+poly_degree-2,nqp));
elseif poly_degree == 1
        basis_functions = struct('N_0',zeros(n+poly_degree,nqp), ...
                                 'N_1',zeros(n+poly_degree-1,nqp));
else
       error('Change polynomial degree!')
end % if

N_final         = zeros(n,nqp);

for qp = 1:nqp
    
    %% constant basis functions 
    p = 0;  
    % determine current number of basis functions n
    n = nkn-p-1;
    % initialization
    N_0 = zeros(n,1);
    % calculation of basis functions
    for i = 1:n
      % check value of xi
      % if xi is the right boundary of the knot vector,
      % then pertubate xi so that boundary value can be determined
        xi_end = 0;
        if xi(qp) == XI(end) 
            xi_2 = xi(qp)-10^(-6); 
            % set N to "1" for specific interval
            if( xi_2 >= XI(i) && xi_2 < XI(i+1) )
                N_0(i,1) = 1;
            end
            xi_end = 1;
        end
        % set N to "1" for specific interval
        if( xi_end ==0 && xi(qp) >= XI(i) && xi(qp) < XI(i+1) )
            N_0(i,1) = 1;
        end
    end
    
    basis_functions.N_0(:,qp) = N_0;
    N = N_0;

    %% linear basis functions   
    if poly_degree > 0
    p = 1;  
    % determine current number of basis functions n
    n = nkn-p-1;
    % initialization
    N_1      = zeros(n,1);
    
      for i = 1:n
        % calculation of the basis functions
        diff1 = XI(i+p) - XI(i);
        diff2 = XI(i+p+1) - XI(i+1);
        % avoid division by zero 
        % source: script "Structure Optimzation"
        if diff1 == 0  
          summand1 = 0;
        else 
          summand1 = (xi(qp) - XI(i))/diff1;
        end
        if diff2 == 0
          summand2 = 0; 
        else
          summand2 = (XI(i+p+1) - xi(qp))/diff2;    
        end
 
        N_1(i,1) = summand1 * N_0(i,1) + summand2 * N_0(i+1,1);
        
      end %for

      basis_functions.N_1(:,qp) = N_1;

      N = N_1;
      
    end % if

    %% quadratic basis functions   
    if poly_degree > 1
    p = 2;  
    % determine current number of basis functions n
    n = nkn-p-1;
    % initialization
    N_2    = zeros(n,1);

      for i = 1:n
        % calculation of the basis functions
        diff1 = XI(i+p) - XI(i);
        diff2 = XI(i+p+1) - XI(i+1);
        % avoid division by zero 
        % source: script "Structure Optimzation"
        if diff1 == 0  
          summand1 = 0;
        else 
          summand1 = (xi(qp) - XI(i))/diff1;
        end
        if diff2 == 0
          summand2 = 0; 
        else
          summand2 = (XI(i+p+1) - xi(qp))/diff2;    
        end

        N_2(i,1) = summand1 * N_1(i,1) + summand2 * N_1(i+1,1);
        
      end %for

      basis_functions.N_2(:,qp) = N_2;

      N = N_2;
    
    end %if

    %% cubic basis functions   
    if poly_degree > 2
    p = 3;  
    % determine number of basis functions n
    n = nkn-p-1;
    % initialization
    N_3    = zeros(n,1);

      for i = 1:n
        % calculation of the basis functions
        diff1 = XI(i+p) - XI(i);
        diff2 = XI(i+p+1) - XI(i+1);
        % avoid division by zero
        % source: script "Structure Optimzation"
        if diff1 == 0  
          summand1 = 0;
        else 
         summand1 = (xi(qp) - XI(i))/diff1;
        end
        if diff2 == 0
          summand2 = 0; 
        else
          summand2 = (XI(i+p+1) - xi(qp))/diff2;    
        end
        N_3(i,1) = summand1 * N_2(i,1) + summand2 * N_2(i+1,1);
        
      end %for

      basis_functions.N_3(:,qp) = N_3;

      N = N_3;
      
    end %if
    
    %% basis functions with p=4 
    if poly_degree > 3
    p = 4;  
    % determine number of basis functions n
    n = nkn-p-1;
    % initialization
    N_4    = zeros(n,1);

      for i = 1:n
        % calculation of the basis functions
        diff1 = XI(i+p) - XI(i);
        diff2 = XI(i+p+1) - XI(i+1);
        % avoid division by zero
        % source: script "Structure Optimzation"
        if diff1 == 0  
          summand1 = 0;
        else 
         summand1 = (xi(qp) - XI(i))/diff1;
        end
        if diff2 == 0
          summand2 = 0; 
        else
          summand2 = (XI(i+p+1) - xi(qp))/diff2;    
        end
        N_4(i,1) = summand1 * N_3(i,1) + summand2 * N_3(i+1,1);
        
      end %for

      basis_functions.N_4(:,qp) = N_4;

      N = N_4;
    
    end %if

    %assembly
    N_final(:,qp)      = N;
 
end % for

% Compute final number of basis functions n
n = nkn-poly_degree-1;

end % function
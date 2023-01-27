function [N,n,dN_dxi] = BSpline_BasisFunctions_for_plot(XI,xi,nkn,poly_degree)
% ---------------------------------------------------------------------
% Subroutine BSpline_BasisFunctions_for_plot.m
% calculation of basis functions at location xi
% and their derivatives w.r.t. parametric coordinate xi
% for a polynomial degree up to 4
%
% Author:           Carina Witt
% Date  :           08.05.2018
%
% Input:    xI              - coordinate
%           XI              - knot vector
%           nkn             - number of knots in the knot vector XI
%           poly_degree     - polynomial degree of the basis functions
%
% Output:   N               - value of basis functions at location xi
%           n               - number of Basis functions
%           dN_dxi          - derivative of basis functions w.r.t. xi
%----------------------------------------------------------------------
    
%%% recursive calculation of the basis functions %%%
%%%  and calculation of their derivative for p<0 %%%
% literature: Cotrell, (2.1 - 2.2) 
% and NURBS book (2.9)

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
  if xi == XI(end) 
    xi_2 = xi-10^(-6); 
    % set N to "1" for specific interval
    if( xi_2 >= XI(i) && xi_2 < XI(i+1) )
        N_0(i,1) = 1;
    end
    xi_end = 1;
  end
  % set N to "1" for specific interval
  if( xi_end ==0 && xi >= XI(i) && xi < XI(i+1) )
    N_0(i,1) = 1;
  end
end
N = N_0;

%% linear basis functions   
if poly_degree > 0
    
    p = 1;  
    % determine current number of basis functions n
    n = nkn-p-1;
    % initialization
    N_1    = zeros(n,1);
    dN_dxi = zeros(n,1);

    for i = 1:n
        % calculation of the basis functions
        diff1 = XI(i+p) - XI(i);
        diff2 = XI(i+p+1) - XI(i+1);
        % avoid division by zero 
        % source: script "Structure Optimzation"
        if diff1 == 0  
          summand1 = 0;
        else 
          summand1 = (xi - XI(i))/diff1;
        end
        if diff2 == 0
          summand2 = 0; 
        else
          summand2 = (XI(i+p+1) - xi)/diff2;    
        end
 
        N_1(i,1) = summand1 * N_0(i,1) + summand2 * N_0(i+1,1);
   
        % calculation of the derivatives w.r.t. xi
        if diff1 == 0  
          summand1 = 0;
        else 
          summand1 = p/diff1;
        end
        if diff2 == 0
          summand2 = 0; 
        else
          summand2 = p/diff2;    
        end

        dN_dxi(i,1) = summand1*N_0(i,1) - summand2*N_0(i+1,1); 
    
    end %for
    
    N = N_1;
     
end % if
  
%% quadratic basis functions   
if poly_degree > 1
    
    p = 2;  
    % determine current number of basis functions n
    n = nkn-p-1;
    % initialization
    N_2    = zeros(n,1);
    dN_dxi = zeros(n,1);

    for i = 1:n
        % calculation of the basis functions
        diff1 = XI(i+p) - XI(i);
        diff2 = XI(i+p+1) - XI(i+1);
        % avoid division by zero 
        % source: script "Structure Optimzation"
        if diff1 == 0  
          summand1 = 0;
        else 
          summand1 = (xi - XI(i))/diff1;
        end
        if diff2 == 0
          summand2 = 0; 
        else
          summand2 = (XI(i+p+1) - xi)/diff2;    
        end

        N_2(i,1) = summand1 * N_1(i,1) + summand2 * N_1(i+1,1);

        % calculation of the derivatives w.r.t. xi
        if diff1 == 0  
          summand1 = 0;
        else 
          summand1 = p/diff1;
        end
        if diff2 == 0
          summand2 = 0; 
        else
          summand2 = p/diff2;    
        end

        dN_dxi(i,1) = summand1*N_1(i,1) - summand2*N_1(i+1,1);
    
    end %for

    N = N_2;
    
end %if
  
%% cubic basis functions   
if poly_degree > 2
    
    p = 3;  
    % determine number of basis functions n
    n = nkn-p-1;
    % initialization
    N_3    = zeros(n,1);
    dN_dxi = zeros(n,1);

      for i = 1:n
        % calculation of the basis functions
        diff1 = XI(i+p) - XI(i);
        diff2 = XI(i+p+1) - XI(i+1);
        % avoid division by zero
        % source: script "Structure Optimzation"
        if diff1 == 0  
          summand1 = 0;
        else 
         summand1 = (xi - XI(i))/diff1;
        end
        if diff2 == 0
          summand2 = 0; 
        else 
          summand2 = (XI(i+p+1) - xi)/diff2;    
        end
        N_3(i,1) = summand1 * N_2(i,1) + summand2 * N_2(i+1,1);

        % calculation of the derivatives w.r.t. xi
        if diff1 == 0  
          summand1 = 0;
        else 
          summand1 = p/diff1;
        end
        if diff2 == 0
          summand2 = 0; 
        else
          summand2 = p/diff2;    
        end

        dN_dxi(i,1) = summand1*N_2(i,1) - summand2*N_2(i+1,1);
      
      end %for

      N = N_3;
end %if

%% basis functions with p=4
if poly_degree > 3
    
    p = 4;  
    % determine number of basis functions n
    n = nkn-p-1;
    % initialization
    N_4    = zeros(n,1);
    dN_dxi = zeros(n,1);

      for i = 1:n
        % calculation of the basis functions
        diff1 = XI(i+p) - XI(i);
        diff2 = XI(i+p+1) - XI(i+1);
        % avoid division by zero
        % source: script "Structure Optimzation"
        if diff1 == 0  
          summand1 = 0;
        else 
         summand1 = (xi - XI(i))/diff1;
        end
        if diff2 == 0
          summand2 = 0; 
        else 
          summand2 = (XI(i+p+1) - xi)/diff2;    
        end
        N_4(i,1) = summand1 * N_3(i,1) + summand2 * N_3(i+1,1);

        % calculation of the derivatives w.r.t. xi
        if diff1 == 0  
          summand1 = 0;
        else 
          summand1 = p/diff1;
        end
        if diff2 == 0
          summand2 = 0; 
        else
          summand2 = p/diff2;    
        end

        dN_dxi(i,1) = summand1*N_3(i,1) - summand2*N_3(i+1,1);
      
      end %for

      N = N_4;
end %if

% Compute final number of basis functions n
n = nkn-poly_degree-1;

end % function
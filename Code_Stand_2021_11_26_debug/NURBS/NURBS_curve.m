function [C] = NURBS_curve(N,n,m,B,w,ndm)
% ---------------------------------------------------------------------
% Subroutine NURBS_Curve.m
% for exact geometry description
% on the basis of B-Spline basis functions and the control polygon
%
% Author:           Carina Witt
% Date  :           15.05.2018
%
% Input:    N               - B-Spline basis functions
%           n,m             - number of basis functions
%           B               - Control points
%           w               - Weights of the control points
%           ndm             - number of dimensions
%
% Output:   C               - NURBS curve
%---------------------------------------------------------------------- 

% literature: Cottrell, (2.28)

% initialization
R = zeros(n,1);
C = zeros(m,ndm);

% loop over all curves
for curve = 1:m
  % Rational basis functions
  sum = 0;
  for j=1:n
    sum = sum + N(j,1)*w(j,1,curve);
  end
  % sum
  for i=1:n
    R(i,1) = N(i,1)*w(i,1,curve)/sum;
  end
  % R
  % NURBS Curve
  for i=1:n
    C(curve,:) = C(curve,:) + R(i,1)*B(i,:,curve);
  end
end
% C

end % function
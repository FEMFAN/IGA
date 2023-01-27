function [A] = NURBS_surface(N,n,M,m,B,w,ndm)
% ---------------------------------------------------------------------
% Subroutine NURBS_Surface.m
% for exact geometry description
% on the basis of B-Spline basis functions and the control polygon
%
% Author:           Carina Witt
% Date  :           15.05.2018
%
% Input:    N,M             - B-Spline basis functions
%           n,m             - number of B-Spline basis functions
%           B               - Control points
%           w               - Weights of the control points
%           ndm             - number of dimensions
%
% Output:   A               - NURBS surface 
%-----------------------------------------------------------------------

% initialization
R = zeros(n,m);
A = zeros(1,ndm);

% Rational basis functions
sum = 0;
for i=1:n
  for j=1:m
    sum = sum + N(i,1)*M(j,1)*w(i,1,j);
  end
end
for i=1:n
  for j=1:m
    R(i,j) = N(i,1)*M(j,1)*w(i,1,j)/sum;
  end
end
% R

% Surface
for i=1:n
  for j=1:m
     A = A + R(i,j)*B(i,:,j);
  end
end
% A

end % function
function [nkn_XI,nekn_XI,XI_elem] = element_extraction(XI)
% ---------------------------------------------------------------------
% Subroutine element_extraction.m
% extract elements as knot spans from the knot vector
%
% Author:           Carina Witt
% Date  :           17.05.2018
%
% Input:    XI                  - knot vector
%
% Output:   nkn_XI              - number of knots in the knot vector
%           nekn_XI             - number of elements (knot spans) in the knot vector
%           XI_elem             - elements of the knot vector
%---------------------------------------------------------------------- 

% number of knots nkn in the knot vector XI
nkn_XI = size(XI,2);

% extract number of knot spans nekn from knot vector
% -> one knot span defines an element
nekn_XI = 0;
for i=2:size(XI,2)
  if XI(i)~=XI(i-1)
    nekn_XI = nekn_XI + 1;
  end
end
%nekn_XI
  
% save knot spans in a vector XI_elem
k=0;
XI_elem = XI(1);
for i=2:size(XI,2)
  if XI(i)~=XI(i-1)
    k=k+1;
    XI_elem(k+1) = XI(i);
  end
end
%XI_elem

end %function
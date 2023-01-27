function [IEN,INC] = connectivity_arrays(XI,ETA,p,q,n,m,ndm,nel,nen,nnp)
% ---------------------------------------------------------------------
% Subroutine connectivity_arrays.m
% determines the connectivity arrays
% INC: global function number + coordinate direction -> NURBS coordinates
% IEN: local function number + element number -> global function number
%
% Author:           Carina Witt
% Date  :           09.07.2018
%
% Input:    XI,ETA              - knot vectors
%           p,q                 - polynomial degrees
%           n,m                 - number of BSpline basis functions
%           ndm                 - number of dimensions
%           nel                 - number of elements
%           nen                 - number of local basis functions
%           nnp                 - number of global basis functions
%
% Output:   IEN, INC            - connectivity arrays
%         
%---------------------------------------------------------------------- 

% literature: Cottrell, Algorithm 7

% initialization of IEN and INC array
INC = zeros(nnp,ndm);
IEN = zeros(nen,nel);

%initialization of local variables
e=0;      %element number
A=0; B=0; %global function numbers
b=0;      %local function number

% construction of connectivity arrays
for j=1:m
  for i=1:n
    
      % update global function number A
      A = A+1;
      
      % assign NURBS coordinates to INC array
      INC(A,1) = i;
      INC(A,2) = j;
      
      if ( i>=(p+1) && j>=(q+1) ) ...
          && ( XI(i)~=XI(i+1) && ETA(j)~=ETA(j+1) ) 
          %07.11.2018: modification for multiple internal knots
          
          % new element
          e=e+1;
          
          for jloc=0:q
            for iloc=0:p
              
                % global function number B
                B = A - jloc*n - iloc;
                % local function number b
                b = jloc*(p+1) + iloc + 1;
                % assign global function number to IEN array
                IEN(b,e) = B;
           
            end %for
          end %for
          
      end %if
      
  end %for
end %for

end % function
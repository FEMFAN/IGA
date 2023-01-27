function [XI_mod,nkn_XI_mod,nekn_XI_mod,XI_elem_mod,B_mod,w_mod,n_mod] = knot_insertion_XI(n,m,ndm,XI,nkn,poly_degree,B,w,new_knot)
% ---------------------------------------------------------------------
% Subroutine knot_insertion_XI.m
% Defines new knot vector and control points
% on the basis of the new knots to be inserted
%
% Author:           Carina Witt
% Date  :           28.05.2018
%                   20.11.2019 modification for multiple interior knots
%                   03.04.2020 modification/correction -> compute on
%                   hyperplane
%
% Input:    n,m             - number of shape functions
%           ndm             - number of dimensions
%           XI              - original knot vector
%           xi              - current coordinate value
%           nkn             - original number of knots
%           poly_degree     - ploynomial degree of the basis functions
%           B               - original control points
%           w               - original weights of the control points
%           new_knot        - knot to be inserted
%
% Output:   XI_mod          - refined knot vector
%           nkn_XI_mod      - new number of knots
%           nekn_XI_mod     - new number of elements (knot spans)
%           XI_elem_mod     - new elements (knot spans)
%           B_mod           - modified control points for new knot vector
%           w_mod           - modified weights
%           n_mod           - new number of basis functions 
%---------------------------------------------------------------------- 

%% Insert new knots into the knot vector
% length of the knot vector is increased by 1

% initialization
XI_mod = XI(1,1);
k = 0;
multiplicity = 0;
% create new knot vector 
for i=2:(nkn-1)
  % copy old knot value to the new knot vector
  XI_mod(1,end+1) = XI(1,i);
  % check position of the new knot
  if ( new_knot >= XI(1,i) && new_knot < XI(1,i+1) )
    % save position of the new knot as k
    k = i;
    % insert new knot into the new knot vector
    XI_mod(1,end+1) = new_knot;
  end 
  
  % check if knot is already present at the prior location the knot vector
  % -> increase of the multiplicity
  if new_knot == XI(1,i)
      multiplicity = multiplicity+1;
  end
  
  % last entry
  if i==(nkn-1)
    % add new knot if it is equal the 2nd last knot
      if new_knot==XI(1,i+1)
        multiplicity = multiplicity+1; % for open knot vector
        XI_mod(1,end+1) = new_knot;
        k = i+1;
      end 
    % add last knot
    XI_mod(1,end+1) = XI(1,nkn);
  end
end % for
%XI_mod

% warning if multiplicity is larger than polynomial degree!
% i.e. the multiplicity before insertion of the new knot is equal to p
if multiplicity >= poly_degree
    fprintf('Warning: The multiplicity of some knots is larger than the polynomial degree!\n');
end

% extract new number of knots and knot spans from knot vector
% via subroutine element_extraction.m
[nkn_XI_mod,nekn_XI_mod,XI_elem_mod] = element_extraction(XI_mod);

%% Compute new control points B and weights w
% continuity of the curve is preserved by computing the control points as below
% source: http://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/NURBS-knot-insert.html
% and: https://pages.mtu.edu/~shene/COURSES/cs3621/LAB/surface/knot-insrt.html

% initialization
n_mod = n;
B_mod = zeros(n_mod,ndm,m);
w_mod = zeros(n_mod,1,m);

% compute weighted control points
Bw_mod = zeros(n_mod,ndm+1,m);
Bw = B(:,:,:).*w(:,1,:);
Bw(:,end+1,:) = w(:,1,:);

% increase number of shape functions by one
% only one knot is inserted at a time!
n_mod=n_mod+1;
for t=1:m
    for i=1:n_mod 
         % special case: new knot is equal to last knot
         if i==n_mod && new_knot==XI(1,i+1) 
           %alpha = 0; 
           Bw_mod(i,:,t) = Bw(i-1,:,t);
           break
         end
         if ( i<=(k-poly_degree) )
           %alpha = 1;
           Bw_mod(i,:,t) = Bw(i,:,t);
         elseif ( i>=(k-poly_degree+1) && i<=(k-multiplicity) )
           alpha = (new_knot-XI(i))/(XI(i+poly_degree)-XI(i));
           Bw_mod(i,:,t) = alpha*Bw(i,:,t)+(1-alpha)*Bw(i-1,:,t);
         elseif ( i>=(k-multiplicity+1) )
           %alpha = 0; 
           Bw_mod(i,:,t) = Bw(i-1,:,t);
         end
    end
end

% update B and w before the next knot is inserted
B_mod = Bw_mod(:,1:ndm,:)./Bw_mod(:,ndm+1,:);
w_mod = Bw_mod(:,ndm+1,:);

end % function

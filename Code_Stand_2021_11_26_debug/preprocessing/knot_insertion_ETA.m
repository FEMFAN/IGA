function [ETA_mod,nkn_ETA_mod,nekn_ETA_mod,ETA_elem_mod,B_mod,w_mod,m_mod] = knot_insertion_ETA(n,m,ndm,ETA,nkn,poly_degree,B,w,new_knot)
% ---------------------------------------------------------------------
% Subroutine knot_insertion_ETA.m
% Defines new knot vector and control points 
% on the basis of the new knots to be inserted in eta-direction
%
% Author:           Carina Witt
% Date  :           18.06.2018
%                   20.11.2019 modification for multiple interior knots
%                   03.04.2020 modification/correction -> compute on
%                   hyperplane
%
% Input:    n,m             - number of basis functions
%           ndm             - number of dimensions
%           ETA             - original knot vector
%           eta             - current coordinate value
%           nkn             - original number of knots
%           poly_degree     - polynomial degree of the basis functions
%           B               - original control points
%           w               - original weights of the control points
%           new_knot        - knot to be inserted
%
% Output:   ETA_mod         - refined knot vector
%           nkn_mod         - new number of knots
%           nekn_mod        - new number of elements (knot spans)
%           ETA_elem_mod    - new elements (knot spans)
%           B_mod           - modified control points for new knot vector
%           w_mod           - modified weights for new knot vector
%           m_mod           - new number of basis functions
%---------------------------------------------------------------------- 

%% Insert new knots into the knot vector
% length of the knot vector is increased by 1 with each inserted knot

% initialization
ETA_mod = ETA(1,1);
k=0;
multiplicity = 0;
% create new knot vector 
for i=2:(nkn-1)
  % copy old knot value to the new knot vector
  ETA_mod(1,end+1) = ETA(1,i);
  % check position of the new knot
  if (new_knot>=ETA(1,i) && new_knot<ETA(1,i+1))
    % save position of the new knot in a vector k
    k = i;
    % insert new knot into the new knot vector
    ETA_mod(1,end+1) = new_knot;
  end 
  
  % check if knot is already present at the prior location the knot vector
  % -> increase of the multiplicity
  if new_knot == ETA(1,i)
      multiplicity = multiplicity+1;
  end
  
  % last entry
  if i==(nkn-1)
    if new_knot==ETA(1,i+1)
      multiplicity = multiplicity+1; % for open knot vector
      ETA_mod(1,end+1) = new_knot;
      k = i+1;
    end 
    ETA_mod(1,end+1) = ETA(1,nkn);
  end
end
%ETA_mod

% warning if multiplicity is larger than polynomial degree!
% i.e. the multiplicity before insertion of the new knot is equal to p
if multiplicity >= poly_degree
    fprintf('Warning: The multiplicity of some knots is larger than the polynomial degree!\n');
end

% extract new number of knots and knot spans from knot vector
% via subroutine element_extraction.m
[nkn_ETA_mod,nekn_ETA_mod,ETA_elem_mod] = element_extraction(ETA_mod);

%% Compute new control points B and weihgts w
% continuity of the curve is preserved by computing the control points as below
% source: http://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/ETARBS-knot-insert.html
% and: https://pages.mtu.edu/~shene/COURSES/cs3621/LAB/surface/knot-insrt.html

% reshape B and w to achieve connectivity in eta-direction
B_re = zeros(m,ndm,n);
w_re = zeros(m,1,n);
for i=1:n
  for j=1:m
  B_re(j,:,i) = B(i,:,j);
  w_re(j,1,i) = w(i,1,j);
  end
end
% compute weighted control points
Bw_re = B_re(:,:,:).*w_re(:,1,:);
Bw_re(:,end+1,:) = w_re(:,1,:);

% initialization
m_mod    = m;
B_re_mod = zeros(m_mod,ndm,n);
w_re_mod = zeros(m_mod,1,n);
Bw_re_mod   = zeros(m_mod,ndm+1,n);

% increase number of shape functions by one
% only one knot is inserted at a time!
m_mod=m_mod+1;
for t=1:n
   for i=1:m_mod 
       % special case: new knot is equal to last knot
       if i==m_mod && new_knot==ETA(1,i+1) 
           %alpha = 0; 
           Bw_re_mod(i,:,t) = Bw_re(i-1,:,t);
           break
       end
       if ( i<=(k-poly_degree) )
         %alpha = 1;
         Bw_re_mod(i,:,t) = Bw_re(i,:,t);
       elseif ( i>=(k-poly_degree+1) && i<=(k-multiplicity) )
         alpha = (new_knot-ETA(i))/(ETA(i+poly_degree)-ETA(i));
         Bw_re_mod(i,:,t) = alpha*Bw_re(i,:,t)+(1-alpha)*Bw_re(i-1,:,t);
       elseif ( i>=(k-multiplicity+1) )
         %alpha = 0; 
         Bw_re_mod(i,:,t) = Bw_re(i-1,:,t);
       end
   end
end
  
% update B and w before the next knot is inserted
B_re = Bw_re_mod(:,1:ndm,:)./Bw_re_mod(:,ndm+1,:);
w_re = Bw_re_mod(:,ndm+1,:);

B_mod = zeros(n,ndm,m_mod);
w_mod = zeros(n,1,m_mod);
for i=1:m_mod
  for j=1:n
  B_mod(j,:,i) = B_re(i,:,j);
  w_mod(j,1,i) = w_re(i,1,j);
  end
end

end % function

function Plot_NURBS_surfaces(XI,ETA,nkn_XI,nkn_ETA,poly_degree,B,w,ndm)
% ---------------------------------------------------------------------
% Subroutine Plot_NURBS_surfaces.m
% creates 2D Plot of the NURBS surface
%
% Author:           Carina Witt
% Date  :           18.05.2018
%
% Input:    XI,ETA          - knot vectors
%           nkn_XI,nkn_ETA  - number of knots in the knot vectors
%           n,m             - number of shape functions
%           poly_degree     - polynomial degree of the basis functions
%           B               - Control points
%           w               - Weights of the control points
%           ndm             - number of dimensions
%
% Output:   Plot of the ETARBS surface
%---------------------------------------------------------------------- 

% determine step size for evaluation
stepsize = 0.05;

% number of evaluation points
evaluation_points = (XI(1,end)-XI(1,1))/stepsize + 1;

% initialization
Surface = zeros(evaluation_points,ndm);
step    = 0;

% evaluation of the surface coordinates at points [xi,eta] 
% with defined step size
for xi=XI(1,1):stepsize:XI(1,end)
  for nu=ETA(1,1):stepsize:ETA(1,end)
   step = step+1;
   
   % Call subfunction to calculate value of the basis functions 
   % at points xi and eta
   [N,n] = BSpline_BasisFunctions_for_plot(XI,xi,nkn_XI,poly_degree);
   [M,m] = BSpline_BasisFunctions_for_plot(ETA,nu,nkn_ETA,poly_degree);
    
   % Call subfunction to calculate value of the surface vector 
   % at points xi,eta
   % format of A: [1,ndm] (row vector)
   [A] = NURBS_surface(N,n,M,m,B,w,ndm);
  
   % save value of the surface vector at current step
   Surface(step,:)=A;
 end
end
%Surface;

% plot surface in physical space
plot(Surface(:,1),Surface(:,2))

% alternative for plot:
%patch(Surface(:,1),Surface(:,2))

end % function
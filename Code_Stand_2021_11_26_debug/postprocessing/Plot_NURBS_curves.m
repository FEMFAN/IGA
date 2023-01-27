function Plot_NURBS_curves(XI,nkn,n,m,poly_degree,B,w,ndm)
% ---------------------------------------------------------------------
% Subroutine Plot_NURBS_curves.m
% creates 2D Plot of the NURBS curves
%
% Author:           Carina Witt
% Date  :           18.05.2018
%
% Input:    XI              - knot vectors
%           nkn             - number of knots in the knot vector XI
%           n,m             - number of shape functions
%           poly_degree     - polynomial degree of the basis functions
%           B               - Control points
%           w               - Weights of the control points
%           ndm             - number of dimensions
%
% Output:   Plot of NURBS curves
%---------------------------------------------------------------------- 

% determine step size for evaluation
stepsize = 0.001;

% number of evaluation points
evaluation_points = (XI(1,end)-XI(1,1))/stepsize + 1;

% initialization
Curve = zeros(evaluation_points,ndm,n);
step  = 0;

% evaluation of the shape functions at points xi with defined step size
for xi=XI(1,1):stepsize:XI(1,end)
 step = step+1;
 
 % call subfunction to calculate value of the basis functions at points xi
 [N,n] = BSpline_BasisFunctions_for_plot(XI,xi,nkn,poly_degree);
  
 % call subfunction to calculate value of the NURBS curve at point xi   
 % format of C: [m,ndm]
 [C] = NURBS_curve(N,n,m,B,w,ndm);
    
 % loop over all m curves
   for k =1:m
     % save value of the coordinates
     Curve(step,:,k)=C(k,:);
   end
end
%Curve;

% plot curves over xi (in parameter space)
color = 'rgbcmykrgbcmyk';
for myi=1:size(Curve,3)
    hold on
    plot(Curve(:,1,myi),Curve(:,2,myi),color(myi))
end
    hold off

end % function
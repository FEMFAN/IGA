function Plot_BasisFunctions(XI,nkn,poly_degree,n)
% ---------------------------------------------------------------------
% Subroutine Plot_BasisFunctions.m
% creates 2D Plot of all n basis functions over the coordinate xi
%
% Author:           Carina Witt
% Date  :           16.05.2018
%
% Input:    XI              - knot vector
%           nkn             - number of knots in the knot vector XI
%           poly_degree     - polynomial degree of the basis functions
%           n               - number of shape functions
%
% Output:   Plot
%---------------------------------------------------------------------- 

% determine step size for evaluation
stepsize = 0.001;

% number of evaluation points
evaluation_points = (XI(1,end)-XI(1,1))/stepsize + 1;

% initialization
Shape = zeros(evaluation_points,2,n);
step  = 0;

% evaluation of the shape functions at points xi with defined step size
for xi=XI(1,1):stepsize:XI(1,end)
    
  step = step+1;
  % save value of xi
  Shape(step,1,:) = xi;
    
  % Call subfunction to calculate value of all basis functions at point xi
  [N,~] = BSpline_BasisFunctions_for_plot(XI,xi,nkn,poly_degree);

  % save calculated values of the n basis functions
  for k =1:n 
    Shape(step,2,k) = N(k,1);
  end
  
end

% plot all n shape functions over the coordinate xi
color = 'rgbcmykrgbcmyk';
for myi=1:size(Shape,3)
    hold on
    plot(Shape(:,1,myi),Shape(:,2,myi),color(myi))
end
    hold off

end % function

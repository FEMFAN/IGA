function Plot_surface_mesh(XI,ETA,XI_elem,ETA_elem,nkn_XI,nkn_ETA,poly_degree,B,w,ndm)
% ---------------------------------------------------------------------
% Subroutine Plot_surface_mesh.m
% creates 2D Plot of the FE mesh
% elements defined by knot spans
%
% Author:           Carina Witt
% Date  :           19.06.2018
%
% Input:    XI,ETA            - knot vectors
%           XI_elem,ETA_elem  - element decomposition of the knot vectors
%           nkn_XI,nkn_ETA    - number of knots in the knot vectors
%           poly_degree       - polynomial degree of the basis functions
%           B                 - Control points
%           w                 - Weights of the control points
%           ndm               - number of dimensions
%
% Output:   Plot of the FE mesh
%---------------------------------------------------------------------- 

% determine step size for evaluation of surface points
stepsize = 0.01;
% number of evaluation points in both directions
evaluation_points_XI = (XI(1,end)-XI(1,1))/stepsize;
evaluation_points_ETA = (ETA(1,end)-ETA(1,1))/stepsize;

%% Evaluation of surface points

%% evaluation of the surface points at points [xi,eta] with defined step size
% xi direction

% initialization
Surface = zeros(1,ndm);
step=0;

myi=1;
while myi<=size(ETA_elem,2)
   eta=ETA_elem(myi);
   
  for xi=XI_elem(1,1):stepsize:XI_elem(1,end)
   step = step+1;
   
   % Call subfunction to calculate value of the basis functions 
   % at points xi and eta
   [N,n] = BSpline_BasisFunctions_for_plot(XI,xi,nkn_XI,poly_degree);
   [M,m] = BSpline_BasisFunctions_for_plot(ETA,eta,nkn_ETA,poly_degree);
    
   % Call subfunction to calculate value of the surface vector at point xi,eta
   % format of A: [1,ndm] (row vector)
   [x_surf] = NURBS_surface(N,n,M,m,B,w,ndm);
  
   % save value of the surface vector
   if step==1
    Surface(1,:)=x_surf;
   else
    Surface(end+1,:)=x_surf;
   end
  end
  myi = myi+1;
end
%end

%% evaluation of the surface points at points [xi,eta] with defined step size
% eta direction

% initialization
Surface_re = zeros(1,ndm);
step=0;

myi=1;
while myi<=size(XI_elem,2)
   xi = XI_elem(myi);
   
   for eta=ETA_elem(1,1):stepsize:ETA_elem(1,end)
       step = step+1;

       % Call subfunction to calculate value of the basis functions 
       % at points xi and eta
       [N,n] = BSpline_BasisFunctions_for_plot(XI,xi,nkn_XI,poly_degree);
       [M,m] = BSpline_BasisFunctions_for_plot(ETA,eta,nkn_ETA,poly_degree);

       % Call subfunction to calculate value of the surface vector at point xi,eta
       % format of A: [1,ndm] (row vector)
       [x_surf] = NURBS_surface(N,n,M,m,B,w,ndm);

       % save value of the surface vector
       if step==1
        Surface_re(1,:)     = x_surf;
       else
        Surface_re(end+1,:) = x_surf;
       end
   end
   
   myi = myi+1;
end

% plot surface points and their connection
% loop is relevant for the right connection of end points
i=1;
while i<= size(Surface,1)-evaluation_points_XI
  plot(Surface(i:(i+evaluation_points_XI),1),Surface(i:(i+evaluation_points_XI),2),'k');
  i = i+evaluation_points_XI+1;
end
hold on
i=1;
while i<= size(Surface_re,1)-evaluation_points_ETA
  plot(Surface_re(i:(i+evaluation_points_ETA),1),Surface_re(i:(i+evaluation_points_ETA),2),'k');
  i = i+evaluation_points_ETA+1;
end

end % function
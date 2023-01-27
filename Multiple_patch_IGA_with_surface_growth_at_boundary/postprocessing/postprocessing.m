function [x,conn] = postprocessing(XI,ETA,XI_elem,ETA_elem,nkn_XI,nekn_XI,nkn_ETA,nekn_ETA,poly_degree,B,w,ndm,nel,n,p)
% ---------------------------------------------------------------------
% Subroutine postprocessing.m
% computation of current physical coordinates and connectivity
% for visualization output (Paraview vtk-file)
%
% Author:           Carina Witt
% Date  :           18.10.2018
%
% Input:    XI,ETA             - knot vectors
%           XI_elem,ETA_elem2 - element decomposition of the knot vectors
%           nkn_XI,nkn_ETA     - number of knots in the knot vectors
%           poly_degree        - polynomial degree of the basis functions
%           B                  - Control points
%           w                  - Weights of the control points
%           ndm                - number of dimensions
%           nel                - number of elements
%           n                  - number of basis funcions in xi-direction
%           p                  - polynomial degree
%
% Output:   x
%           conn
%           Plot of the FE mesh
%---------------------------------------------------------------------- 

% initialization
nodecoordstruct = struct('X',zeros(4,ndm)); %for 4-noded elements
node_coords     = repmat(nodecoordstruct,nel,1);
nodestruct      = struct('X',zeros(1,ndm)); %for 4-noded elements
node            = repmat(nodestruct,nel*4,1);

step_XI  = XI_elem(2)-XI_elem(1);
step_ETA = ETA_elem(2)-ETA_elem(1);

np = 0;
step = 0;      
myi=1;
while myi<=size(ETA_elem,2)
   eta = ETA_elem(myi);
   myj = 1;

   while myj<=size(XI_elem,2)
   xi = XI_elem(myj);

       np=np+1;
       step = step+1;

       % Call subfunction to calculate value of the basis functions 
       % at points xi and eta
       [N,n] = BSpline_BasisFunctions_for_plot(XI,xi,nkn_XI,poly_degree);
       [M,m] = BSpline_BasisFunctions_for_plot(ETA,eta,nkn_ETA,poly_degree);

       % Call subfunction to calculate value of the surface vector at point xi,eta
       % format of A: [1,ndm] (row vector)
       [x_surf]  = NURBS_surface(N,n,M,m,B,w,ndm);
       x(step,:) = x_surf;

       myj = myj+1;

  end

   myi = myi+1;

end

% connectivity list
conn = zeros(nel,4);
i=0;
flag=0;
for e=1:nel
    i=i+1;
    if flag==1
        i=i+1;
        flag = 0;
    end
    a1 = i;
    a2 = i+1;
    c1 = 2*(nekn_XI+1)-nekn_XI+i-1;
    c2 = 2*(nekn_XI+1)-nekn_XI+i;  
    
    if mod(c2,(n-p+1))==0
        flag = 1;
    end

    conn(e,:) = [a1 a2 c2 c1];
end


end % function
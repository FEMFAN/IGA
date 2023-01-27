function Plot_controlpolygon(B,w,n,m,ndm)
% ---------------------------------------------------------------------
% Subroutine Plot_controlpolygon.m
% creates 2D Plot of the control polygon
%
% Author:           Carina Witt
% Date  :           14.06.2018
%
% Input:    B              - Control Points
%           w              - weights of the control points
%           n,m            - number of basis functions
%           ndm            - number of dimensions (here: ndm=2)
%
% Output:   Plot of control polygon
%---------------------------------------------------------------------- 

% reshape B for plot to achieve connectivity in nu-direction
B_re = zeros(m,ndm,n);
for i=1:n
  for j=1:m
    B_re(j,:,i) = B(i,:,j);
  end
end

% 2D plot of control polygon

% plot control points and their connection in xi-direction
for i=1:(m)
    plot(B(:,1,i),B(:,2,i),'-ok');
    hold on
end
% plot connections of the control points in nu-direction
for j=1:n
    plot(B_re(:,1,j),B_re(:,2,j),'k');
    hold on
end

hold off

end % function
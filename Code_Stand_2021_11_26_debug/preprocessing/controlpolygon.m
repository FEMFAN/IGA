function [XI,ETA,nkn_XI,nekn_XI,XI_elem,nkn_ETA,nekn_ETA,ETA_elem,KP,w8,B,w,n,m] = controlpolygon(n,m,p,q,XI,ETA,nkn_XI,nekn_XI,XI_elem,nkn_ETA,nekn_ETA,ETA_elem,ndm,KP,w8,knot_removal_flag,remove_knots_XI,h_refinement_flag,new_knots_XI,new_knots_ETA)
% ---------------------------------------------------------------------
% Subroutine controlpolygon.m
% extract elements as knot spans from the knot vector
% and constructs and plots the control polygon
% reshapes control points and weights
%
% Author:           Carina Witt
% Date  :           17.05.2018
%                   13.03.2020 modification for knot removement
%
% Input:    ...
%
% Output:   XI,ETA                             - extended knot vectors
%           nkn_XI,nekn_XI,nkn_ETA,nekn_ETA    - new element information
%           XI_elem,ETA_elem                   - element information
%           KP,w8                              - refined control polygon
%           B,w                                - control points and weights
%           n,m                                - new numbers of basis funcions
%       
%---------------------------------------------------------------------- 

% reshape list of control points and weights
% extract control points and weights for each individual curve
B = zeros(n,ndm,m);
w = zeros(n,1,m);
for i = 0:(m-1)
  B(:,:,i+1) = KP((i*n+1):((i+1)*n),:);
  w(:,1,i+1) = w8((i*n+1):((i+1)*n),1);
end

% plot control polygon
figure(001)
hold on
Plot_controlpolygon(B,w,n,m,ndm)
axis('equal')
hold off

%% mesh modifications


% if desired: knot removal
% 0 - no removal
% 1 - knot removal active for knot vector XI
% 2 - knot removal active for knot vector ETA
% 3 - knot removal active for both knot vectors

if knot_removal_flag ~= 0

    if knot_removal_flag == 1 || knot_removal_flag == 3
      % overwrite knot vector, shape functions and control points
      [XI,nkn_XI,nekn_XI,XI_elem,B,w,n] = knot_removal_XI(m,ndm,XI,p,B,w,remove_knots_XI);
    end %if

    if knot_removal_flag ==2 || knot_removal_flag == 3
      fprintf('knot removal for ETA not implemented yet!')
    end %if

    % reshape B into KP 
    KP = zeros(n*m,ndm);
    w8 = zeros(n*m,1);
    for i = 0:(m-1)
        KP((i*n+1):((i+1)*n),:) = B(:,:,i+1);
        w8((i*n+1):((i+1)*n),1) = w(:,1,i+1);
    end
  
  
  % plot reduced control polygon 
  figure(002)
  hold on
  grid on
  Plot_controlpolygon(B,w,n,m,ndm)
  axis('equal')
  hold off
  
end %if knot_removal_flag~=0


% if desired: knot insertion (h-refinement)
% 0 - no refinement
% 1 - knot insertion active for knot vector XI
% 2 - knot insertion active for knot vector ETA
% 3 - knot insertion active for both knot vectors

if h_refinement_flag ~= 0

    if h_refinement_flag == 1 || h_refinement_flag == 3
      % loop over all knots to be inserted
      for newknot=1:size(new_knots_XI,2)
          new_knot_XI = new_knots_XI(1,newknot);     
          % overwrite knot vector, shape functions and control points
          [XI,nkn_XI,nekn_XI,XI_elem,B,w,n] = knot_insertion_XI(n,m,ndm,XI,nkn_XI,p,B,w,new_knot_XI);
      end
    end %if

    if h_refinement_flag ==2 || h_refinement_flag == 3
      for newknot=1:size(new_knots_ETA,2)    
          new_knot_ETA = new_knots_ETA(1,newknot);   
          [ETA,nkn_ETA,nekn_ETA,ETA_elem,B,w,m] = knot_insertion_ETA(n,m,ndm,ETA,nkn_ETA,q,B,w,new_knot_ETA);
      end
    end %if

    % reshape B into KP 
    KP = zeros(n*m,ndm);
    w8 = zeros(n*m,1);
    for i = 0:(m-1)
        KP((i*n+1):((i+1)*n),:) = B(:,:,i+1);
        w8((i*n+1):((i+1)*n),1) = w(:,1,i+1);
    end
  
  % plot refined control polygon 
  figure(003)
  hold on
  grid on
  Plot_controlpolygon(B,w,n,m,ndm)
  axis('equal')
  hold off
  
end %if h_refinement_flag~=0

end %function

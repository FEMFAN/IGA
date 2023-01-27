function [ID,drclDofs,freeDofs,neumDofs,fpre,u_d] = ID_array(ndm,nnp,INC,u_pre,f_pre,loadcurve,time)
% ---------------------------------------------------------------------
% Subroutine ID_array.m
% ID: global function number + coordinate direction -> global equation number
% determines also the free, Dirichlet and Neumann dofs
% and the vector of prescribed forces
%
% Author:           Carina Witt
% Date  :           31.07.2018
%
% Input:      ndm               - number of dimensions
%             nnp               - number of global basis functions
%             INC               - INC connectivity array
%             u_pre,f_pre       - prescribed displacement and force
%             loadcurve, time   - loadcurve and current time
%
% Output:     ID                - connectivity array
%             drclDofs          - dofs on Dirichlet boundary 
%             freeDofs          - dofs not on Dirichlet boundary
%             neumDofs          - dofs on Neumann boundary
%             fpre              - prescribed force vector
%             u_d               - prescribed displacement values
%
%---------------------------------------------------------------------- 

% initialization
ID       = zeros(ndm,nnp);
drclDofs = [];
freeDofs = [];
neumDofs = [];
fpre     = zeros(nnp*ndm,1);

% set counters to zero
counter1     = 0;
counter2     = 0;
counter3     = 0;
counter4     = 0;
glob_fun_dim = 0;

% construction of connectivity array and determination of dofs
for glob_fun = 1:nnp
    local_fun_xi  = INC(glob_fun,1);
    local_fun_eta = INC(glob_fun,2);
    
  for dim = 1:ndm
       glob_fun_dim = glob_fun_dim+1;
       flag=0;
       if size(u_pre,1) == 0
           drclDofs = [];
            u_d = [];
       else
       for i=1:size(u_pre,1)
           
           node_XI  = u_pre(i,1);
           node_ETA = u_pre(i,2);
           ldof     = u_pre(i,3);
           loadid   = u_pre(i,4);
           scale    = u_pre(i,5);
           
           % u_pre: [INC direction (xi or nu), INC value (i.e. local basis fct number for the specific direction, ...
           % displacement direction, displacement value]
           if    ( ldof==dim  && local_fun_xi == node_XI ...
              && (node_ETA==0 || local_fun_eta == node_ETA) ) || ...
                 ( ldof==dim  && local_fun_eta == node_ETA ...
              && (node_XI ==0 || local_fun_xi == node_XI) )   || ...
                 ( ldof==dim && node_XI==0 && node_ETA==0 ) %all basis functions %22.10.18
               %Dirichlet boundary (prescribed displacement)
               counter2 = counter2+1;
               % set ID to zero on Dirichlet boundary
               ID(dim,glob_fun)   = 0;
               % update Dirichlet dofs
               drclDofs(counter2) = glob_fun_dim;
               % determine vector of prescribed displacements
               u_d(counter2) = scale * interp1(loadcurve(loadid).time, loadcurve(loadid).value, time);
               flag =1;
           end
           
       end
       end
       
       if flag==0 % not on Dirichlet boundary
           counter1 = counter1+1;
           counter3 = counter3+1;
           % update equation numbering in the ID array
           ID(dim,glob_fun)   = counter1;
           % update free dofs
           freeDofs(counter3) = glob_fun_dim;
       end
       
       % build vector of prescribed forces
       for j=1:size(f_pre,1)
           
           node_XI  = f_pre(j,1);
           node_ETA = f_pre(j,2);
           ldof     = f_pre(j,3);
           loadid   = f_pre(j,4); 
           scale    = f_pre(j,5);
           
           if    ( ldof==dim  && local_fun_xi == node_XI ...
              && (node_ETA==0 || local_fun_eta == node_ETA) ) || ...
                 ( ldof==dim  && local_fun_eta == node_ETA ...
              && (node_XI ==0 || local_fun_xi == node_XI) )  || ...
                 ( ldof==dim  && node_XI==0 && node_ETA==0 ) %all basis functions %22.10.18
              % Neumann boundary (prescribed force)
              counter4 = counter4+1;
              % update Neumann dofs
              neumDofs(counter4) = glob_fun_dim;
              % determine vector of prescribed forces
              fpre(glob_fun_dim) = scale * interp1(loadcurve(loadid).time, loadcurve(loadid).value, time);
           end
           
       end
       
  end
  
end
%ID

% transpose u_d
if size(u_d,1)==0
else
u_d = u_d';
end
end % function
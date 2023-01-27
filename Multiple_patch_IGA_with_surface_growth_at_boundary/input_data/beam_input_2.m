function [fid,u_pre,f_pre,f_b,ttime,dt,loadcurve,D1_flag] = beam_input_2(n,m,p,k,step)

% choose boundary problem
load_flag = 1; %=1: tension test; %=2: shear test
drive_flag = 1; %=1: displacement_driven test; %=2: force_driven test
D1_flag = 0;    %=0: 2 dimensional problem; %=1: one dimensional problem
if load_flag ==1
    fid = 'truss_tipp_growth_BCTipp';
    ldir = 1;
    dvalue = .1;
    fvalue = 600;
elseif load_flag ==2
    fid = 'beam_shear';
    ldir = 2;
    dvalue = 2;
    fvalue = 7000;
end

ttime = 1; % total time of simulation
dt    = .125; % time step

%growth (here solely growth in y direction)
%first row x coordinates of the growth starting point
%second row y coordinates of the growth starting point
%third row x coordinates of the growth ending point
%fourth row y coordinates of the growth ending point
growth = [12    18;
           1    1;
           12   18;
           2    2];

%loadcurves
loadcurve(1).time  =  [0.0   ttime];
loadcurve(1).value = [0.0   1];

% prescribed boundary conditions
% format:
%[1: INC value (i.e. local basis fct number) for the xi direction (if no restriction: =0),...
% 2: INC value (i.e. local basis fct number) for the eta direction (if no restriction: =0)
% 3: displacement/force direction (1=x, 2=y), 4: loadID, 5: displacement/force value]

if load_flag ==1 %tension
    
    if drive_flag ==1
        if k ==1
            u_pre =  [1,1,2,1,0;
                1,0,1,1,0];
            
            
            
            f_pre = [0,0,0,1,0];
            
        elseif k == 2
            
            u_pre = [];
            
            f_pre = [0,0,0,1,0];
        elseif k==3
            u_pre = [n-(8-step),0,ldir,1,dvalue];
            
            f_pre = [0,0,0,1,0];
        end
        
    elseif drive_flag ==2
        if k == 1
            
            u_pre = [1,0,1,1,0;
                1,1,2,1,0];
            f_pre = [0,0,0,1,0];
        elseif k == 2
            
            u_pre = [];
            f_pre = [n,1,ldir,1,fvalue];
            if m==2
                f_pre = [n,1,ldir,1,fvalue/m; %for m=2
                    n,2,ldir,1,fvalue/m;];
            elseif m==3
                f_pre = [n,1,ldir,1,fvalue/m; %for m=3
                    n,2,ldir,1,fvalue/m;
                    n,3,ldir,1,fvalue/m];
            elseif m==4
                f_pre = [n,1,ldir,1,fvalue/(m+2); %for m=4
                    n,2,ldir,1,2*fvalue/(m+2);
                    n,3,ldir,1,2*fvalue/(m+2);
                    n,4,ldir,1,fvalue/(m+2)];
            elseif m>4 && p==1
                for i=1:m
                    if i==1 || i==m
                        j=1;
                    else
                        j=2;
                    end
                    f_pre(i,:) = [n,i,ldir,1,j*fvalue/(2+(m-2)*2)];
                end
            elseif m>4 && p==2
                for i=1:m
                    if i==1 || i==m
                        j=1;
                    elseif i==2 || i==(m-1)
                        j=2;
                    else
                        j=3;
                    end
                    f_pre(i,:) = [n,i,ldir,1,j*fvalue/(6+(m-4)*3)];
                end
            else
                error('Check m for Neumann boundary conditions!');
            end
            f_pre = [n,1,ldir,1,fvalue];
        elseif k ==3
            u_pre = [];
            
            if m==2
                f_pre = [n,1,ldir,1,fvalue/m; %for m=2
                    n,2,ldir,1,fvalue/m;];
            elseif m==3
                f_pre = [n,1,ldir,1,fvalue/m; %for m=3
                    n,2,ldir,1,fvalue/m;
                    n,3,ldir,1,fvalue/m];
            elseif m==4
                f_pre = [n,1,ldir,1,fvalue/(m+2); %for m=4
                    n,2,ldir,1,2*fvalue/(m+2);
                    n,3,ldir,1,2*fvalue/(m+2);
                    n,4,ldir,1,fvalue/(m+2)];
            elseif m>4 && p==1
                for i=1:m
                    if i==1 || i==m
                        j=1;
                    else
                        j=2;
                    end
                    f_pre(i,:) = [n,i,ldir,1,j*fvalue/(2+(m-2)*2)];
                end
            elseif m>4 && p==2
                for i=1:m
                    if i==1 || i==m
                        j=1;
                    elseif i==2 || i==(m-1)
                        j=2;
                    else
                        j=3;
                    end
                    f_pre(i,:) = [n,i,ldir,1,j*fvalue/(6+(m-4)*3)];
                end
            else
                error('Check m for Neumann boundary conditions!');
            end
            %f_pre = [n,1,ldir,1,fvalue];
        end
        
    end
    
elseif load_flag ==2    % shear
    
    if drive_flag ==1
        
        u_pre = [1,0,1,1,0;
            1,0,2,1,0;
            n,0,1,1,0;
            n,0,2,1,dvalue];
        %       for i=2:(n-1) %set x-displacement of middle nodes to zero for homogeneous shear
        %           u_pre(end+1,:) = [i,0,1,1,0];
        %       end
        
        f_pre = [0,0,0,1,0];
        
    elseif drive_flag ==2
        
        u_pre = [1,0,1,1,0;
            1,0,2,1,0;
            n,0,1,1,0];
        for i=2:(n-1) %set x-displacement of middle nodes to zero for homogenous shear
            u_pre(end+1,:) = [i,0,1,1,0];
        end
        
        if m==2
            f_pre = [n,1,ldir,1,fvalue/m; %for m=2
                n,2,ldir,1,fvalue/m;];
        elseif m==3
            f_pre = [n,1,ldir,1,fvalue/m; %for m=3
                n,2,ldir,1,fvalue/m;
                n,3,ldir,1,fvalue/m];
        elseif m==4
            f_pre = [n,1,ldir,1,fvalue/(m+2); %for m=4
                n,2,ldir,1,2*fvalue/(m+2);
                n,3,ldir,1,2*fvalue/(m+2);
                n,4,ldir,1,fvalue/(m+2)];
        elseif m>4
            for i=1:m
                if i==1 || i==m
                    j=1;
                elseif i==2 || i==(m-1)
                    j=2;
                else
                    j=3; %=p+1
                end
                f_pre(i,:) = [n,i,ldir,1,j*fvalue/(6+(m-4)*3)];
            end
        else
            error('Check m for Neumann boundary conditions!');
        end
        
    end
    
end

% body force aka gravity
% value in x direction, value in y direction
f_b = [0 0];

end

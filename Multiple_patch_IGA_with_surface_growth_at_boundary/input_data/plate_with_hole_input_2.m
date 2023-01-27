function [fid,u_pre,f_pre,f_b,ttime,dt,loadcurve,D1_flag] = plate_with_hole_input_2(n,m,p,k)

% choose boundary problem
drive_flag = 1; %=1: displacement_driven test; %=2: force_driven test
D1_flag = 0;    %=0: 2 dimensional problem; %=1: one dimensional problem

fid = 'plate_with_hole_C1_fine';
ldir = 1;
dvalue = .1;
fvalue = 25000;


ttime = 1; % total time of simulation
dt    = 1; % time step

%loadcurves
loadcurve(1).time  =  [0.0   ttime];
loadcurve(1).value = [0.0   1];

% prescribed boundary conditions
% format:
%[1: INC value (i.e. local basis fct number) for the xi direction (if no restriction: =0),...
% 2: INC value (i.e. local basis fct number) for the eta direction (if no restriction: =0)
% 3: displacement/force direction (1=x, 2=y), 4: loadID, 5: displacement/force value]

if drive_flag ==1
    if k ==1
        u_pre =  [1,0,2,1,0];
        u_pre_1 = zeros(floor(n/2),5);
        for i = 1:floor(n/2)
            u_pre_1(i,:) = [i,m,1,1,-dvalue];
        end
        u_pre = [u_pre ; u_pre_1];
        
        f_pre = [0,0,0,1,0];
        
    elseif k == 2
        u_pre =  [n,0,2,1,0];
        u_pre_1 = zeros(floor(n/2),5);
        for i = 1:floor(n/2)
            u_pre_1(i,:) = [i+floor(n/2),m,1,1,dvalue];
        end
        u_pre = [u_pre ; u_pre_1];
        
        f_pre = [0,0,0,1,0];
        
    elseif drive_flag ==2
        if k == 1
            
            u_pre = [1,0,2,1,0];
            
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
                    f_pre(i,:) = [n,i,ldir,1,-j*fvalue/(6+(m-4)*3)];
                end
            else
                error('Check m for Neumann boundary conditions!');
            end
        elseif k == 2
            
            u_pre = [n,0,2,1,0];
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
                for i=1:floor(n/2)
                    if i==1 || i==m
                        j=1;
                    elseif i==2 || i==(m-1)
                        j=2;
                    else
                        j=3;
                    end
                    f_pre(i,:) = [i+floor(n/2),m,ldir,1,j*fvalue/(6+(m-4)*3)];
                end
            else
                error('Check m for Neumann boundary conditions!');
            end
        end
        
        
    end
end
% body force aka gravity
% value in x direction, value in y direction
f_b = [0 0];

end

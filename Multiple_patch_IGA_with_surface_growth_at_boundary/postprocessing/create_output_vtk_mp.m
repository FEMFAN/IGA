function create_output_vtk_mp(filename, np, step, time, ndf, nnp, nel, nqp, X, u, conn, celldata,patches)

%% compute mean stresses and strains on each element
tdm = 2; % tensor dimension
elem=[];
for pl=1:np
    elemstruct       = struct('cn',zeros(1,4),'eps_mean',zeros(tdm*tdm,1),'sig_mean',zeros(tdm*tdm,1), ...
        'kap_mean',zeros(tdm*tdm*tdm,1),'tau_mean',zeros(tdm*tdm*tdm,1));
    elem             = repmat(elemstruct,length(elem)+patches(pl).nel*4,1);
    
end
A = 0;
B = 0;
for pl=1:np
    for e=1:patches(pl).nel
        % save connectivity for element e
        elem(e+A).cn = patches(pl).conn(e,:)+B;
        % calculate mean stresses and strains for element e
        % from the respective gauss point values
        eps_mean = zeros(tdm*tdm,1);
        sig_mean = zeros(tdm*tdm,1);
        
        % reference volume of element e
        V_elem = 0;
        for j=1:nqp
            V_elem = V_elem + patches(pl).celldata(e).dV(j,1);
        end
        
        for i = 1:tdm*tdm
            for j=1:nqp
                eps_mean(i,1) = eps_mean(i,1) + patches(pl).celldata(e).epsilon(i,j) * patches(pl).celldata(e).dV(j,1);
                sig_mean(i,1) = sig_mean(i,1) + patches(pl).celldata(e).sigma(i,j) * patches(pl).celldata(e).dV(j,1);
            end
            eps_mean(i,1) = 1/V_elem * eps_mean(i);
            sig_mean(i,1) = 1/V_elem * sig_mean(i);
            
            % save values
            elem(e+A).eps_mean(i) = eps_mean(i);
            elem(e+A).sig_mean(i) = sig_mean(i);
        end
    end
    B = elem(e+A).cn(end-1);
    A = A + patches(pl).nel;
    
end

%% write vtk-file
% similar to abraxas output file

% filename
vtk_filename = [filename,'_0',num2str(step,'%3d'), '.vtk'];
% write file as text file
fid = fopen(vtk_filename,'Wt');
[fid, message] = fopen(vtk_filename,'Wt');
if fid < 0
    error('Failed to open myfile because: %s', message);
end

% "HEADER" of vtk-file
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, '%s\n', vtk_filename);

fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

% use field data to store time and cycle information
fprintf(fid,'FIELD FieldData 2\n');
fprintf(fid,'TIME 1 1 double\n');
fprintf(fid,'%e\n', time);
fprintf(fid,'CYCLE 1 1 int\n');
fprintf(fid,'%d\n', step);

% "POINTS": node points
X = [];
A=0;
for pl = 1:np
    X=[X; patches(pl).x_initial];
    A = A +patches(pl).nnp_new;
end
fprintf(fid, 'POINTS %d double\n', A);
A = 0;
for pl = 1:np
    for i = 1:patches(pl).nnp_new
        coord        = zeros(1,3);
        coord(1:ndf) = X(i+A,:);
        fprintf(fid, '%24.12e %24.12e %24.12e\n', coord);
    end
    A = A + i;
end

% "CELLS": write element connectivity
size = 0;
A = 0;
for pl = 1:np
    for i = 1:patches(pl).nel
        nen = length(elem(i+A).cn);
        size = size + nen + 1; %+1 f�r zusatzinformation ueber elementtyp
    end
    A = A + patches(pl).nel;
end
fprintf(fid, 'CELLS %d %d\n', A, size);
A=0;
for pl = 1:np
    for i = 1:patches(pl).nel
        nen = length(elem(i+A).cn);
        fprintf(fid, '%d ',nen);  %here: 4
        for j = 1:nen
            % vtk uses C arrays
            % numbering in C starts with 0 not with 1
            fprintf(fid, ' %8d', elem(i+A).cn(j)-1);
        end
        fprintf(fid, '\n');
    end
    A = A +patches(pl).nel;
end

% "CELL_TYPES": write type of each element (BAR2, QUAD4, etc.)
% at the moment only quad4 is used -> vtkElType = 9 see manual
A=0;
for pl=1:np
    A= A +patches(pl).nel;
end
fprintf(fid, 'CELL_TYPES %d\n', A);
A=0;
for pl =1:np
    for i = 1+A:patches(pl).nel+A %+A rausgenommen
        fprintf(fid, '%d\n', 9); % TODO: variable element type
    end
    A= A +patches(pl).nel;
end

% "POINT_DATA": write SCALARS, VECTORS, TENSORS (in that order!)
A=0;
for pl = 1:np
A = A + patches(pl).nnp_new;
end
fprintf(fid, 'POINT_DATA %d\n', A);
% "SCALARS": no SCALARS implemented

% "VECTORS": write displacements
fprintf(fid, 'VECTORS DSPL double\n');
u = [];
for pl = 1:np
    u = [u; patches(pl).u_new];
end

A=0;
for pl=1:np
    for i = 1:patches(pl).nnp_new
        uline        = zeros(1,3);
        uline(1:ndf) = u(i+A,:);
        fprintf(fid, '%24.12e %24.12e %24.12e\n', uline);
    end
    A = A + patches(pl).nnp_new;
end

% "CELL_DATA": write SCALARS, VECTORS, TENSORS (in that order!!)
A=0;
for pl=1:np
    A= A +patches(pl).nel;
end
fprintf(fid, 'CELL_DATA %d\n', A);

% write strains
for j=1:(tdm*tdm)
    if j==1
        varName = 'eps11';
    elseif j==2
        varName = 'eps12';
    elseif j==3
        varName = 'eps21';
    elseif j==4
        varName = 'eps22';
    elseif j==5
        varName = 'eps13';
    elseif j==6
        varName = 'eps23';
    elseif j==7
        varName = 'eps31';
    elseif j==8
        varName = 'eps32';
    elseif j==9
        varName = 'eps33';
    end
    fprintf(fid, 'SCALARS %s double\n', varName);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    A = 0;
    for pl = 1:np
        for i=1:patches(pl).nel
            eps = elem(i+A).eps_mean(j);
            fprintf(fid, '%24.12e\n', eps);
        end
    A = A + patches(pl).nel;
    end
end

% write stresses
for j=1:(tdm*tdm)
    if j==1
        varName = 'sig11';
    elseif j==2
        varName = 'sig12';
    elseif j==3
        varName = 'sig21';
    elseif j==4
        varName = 'sig22';
    elseif j==5
        varName = 'sig13';
    elseif j==6
        varName = 'sig23';
    elseif j==7
        varName = 'sig31';
    elseif j==8
        varName = 'sig32';
    elseif j==9
        varName = 'sig33';
    end
    fprintf(fid, 'SCALARS %s double\n', varName);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    A = 0;
    for pl = 1:np
        for i=1:patches(pl).nel
            sig = elem(i+A).sig_mean(j);
            fprintf(fid, '%24.12e\n', sig);
        end
        A = A + patches(pl).nel;
    end
end

fclose(fid);
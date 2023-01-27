function create_output_vtk(filename, pn, step, time, ndf, nnp, nel, nqp, X, u, conn, celldata)

%% compute mean stresses and strains on each element
tdm = 2; % tensor dimension
elemstruct       = struct('cn',zeros(1,4),'eps_mean',zeros(tdm*tdm,1),'sig_mean',zeros(tdm*tdm,1), ...
                                     'kap_mean',zeros(tdm*tdm*tdm,1),'tau_mean',zeros(tdm*tdm*tdm,1));
elem             = repmat(elemstruct,nel,1);
for e=1:nel
    % save connectivity for element e
    elem(e).cn = conn(e,:);
    % calculate mean stresses and strains for element e
    % from the respective gauss point values
    eps_mean = zeros(tdm*tdm,1);
    sig_mean = zeros(tdm*tdm,1);
    
    % reference volume of element e
    V_elem = 0;
    for j=1:nqp
        V_elem = V_elem + celldata(e).dV(j,1);       
    end

    for i = 1:tdm*tdm
        for j=1:nqp
            eps_mean(i,1) = eps_mean(i,1) + celldata(e).epsilon(i,j) * celldata(e).dV(j,1);
            sig_mean(i,1) = sig_mean(i,1) + celldata(e).sigma(i,j) * celldata(e).dV(j,1);
        end
    eps_mean(i,1) = 1/V_elem * eps_mean(i);
    sig_mean(i,1) = 1/V_elem * sig_mean(i);
    
    % save values
    elem(e).eps_mean(i) = eps_mean(i);
    elem(e).sig_mean(i) = sig_mean(i);
    end

end

%% write vtk-file
% similar to abraxas output file

% filename
vtk_filename = [filename,num2str(pn,'%3d'),'_0',num2str(step,'%3d'), '.vtk'];
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
fprintf(fid, 'POINTS %d double\n', nnp);
for i = 1:nnp
    coord        = zeros(1,3);
    coord(1:ndf) = X(i,:);
    fprintf(fid, '%24.12e %24.12e %24.12e\n', coord);
end

% "CELLS": write element connectivity
size = 0;
for i = 1:nel
    nen = length(elem(i).cn);
    size = size + nen + 1; %+1 f�r zusatzinformation ueber elementtyp
end

fprintf(fid, 'CELLS %d %d\n', nel, size); 
for i = 1:nel
    nen = length(elem(i).cn);
    fprintf(fid, '%d ',nen);  %here: 4
    for j = 1:nen
        % vtk uses C arrays
        % numbering in C starts with 0 not with 1
        fprintf(fid, ' %8d', elem(i).cn(j)-1);
    end
    fprintf(fid, '\n');
end

% "CELL_TYPES": write type of each element (BAR2, QUAD4, etc.)
% at the moment only quad4 is used -> vtkElType = 9 see manual
fprintf(fid, 'CELL_TYPES %d\n', nel);
for i = 1:nel
    fprintf(fid, '%d\n', 9); % TODO: variable element type
end

% "POINT_DATA": write SCALARS, VECTORS, TENSORS (in that order!)
fprintf(fid, 'POINT_DATA %d\n', nnp);

% "SCALARS": no SCALARS implemented

% "VECTORS": write displacements
fprintf(fid, 'VECTORS DSPL double\n');
for i = 1:nnp
    uline        = zeros(1,3);
    uline(1:ndf) = u(i,:);
    fprintf(fid, '%24.12e %24.12e %24.12e\n', uline);
end

% "CELL_DATA": write SCALARS, VECTORS, TENSORS (in that order!!)
fprintf(fid, 'CELL_DATA %d\n', nel);

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
    for i=1:nel
        eps = elem(i).eps_mean(j);
        fprintf(fid, '%24.12e\n', eps);
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
    for i=1:nel
        sig = elem(i).sig_mean(j);
        fprintf(fid, '%24.12e\n', sig);
    end
end

fclose(fid);
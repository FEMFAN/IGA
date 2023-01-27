% Input.m for 2D truss  with no constant cross section
close all;
%clear all;
clc;
np = 2; %number of patches

%interface matrix 1: number of the interface 2%3: number of the meeting
%patches 4: master patch of the interface 5 Interface position inside the
%masters patch interface matrix 6:first coloum of the interface dofs inside the global stiffnes
%matrix 7: last coloum of the interface dofs inside the global stiffnes matrix
IM = [1,1,2,2,1,0,0];

%Patch matrix 1: patch number 2:first coloum of the not in the interface dofs inside the global stiffnes matrix
%3: last coloum of the not in the interface dofs inside the global stiffnes matrix
PM =[1,0,0;
    2,0,0];
%Close to the interface matrix 1: interface number 2:first coloum of the next to the interface dofs inside the global stiffnes matrix
%3: last coloum of the next to the interface dofs inside the global stiffnes matrix
CM = [1,0,0];


%initialize array for patch information
patchstruc = struct('KP',[],'w8',[],'Int',[],'ndm',[],'XI',[],'ETA',[],'nkn_XI',[],'nekn_XI',[],'XI_elem',[],'nkn_ETA',[],'nekn_ETA',[],'ETA_elem',[],'poly_degree',[],'p',[],'q',[],'n',[],'m',[],...
    'B',[],'w',[],'h_refinement_flag',[],'new_knots_XI',[],'new_knots_ETA',[],'nel',[],'nnp',[],'nen',[],'nqp',[],'ndf',[],'tdm',[],'u',[],'K',[],...
    'fint',[],'fext',[],'fvol',[],'frea',[],'u_pre',[],'f_pre',[],'f_b',[],'IEN',[],'INC',[],'drclDofs',[],'freeDofs',[],'fpre',[],...
    'u_d',[],'cellstruct',[],'celldata',[],'u_reshaped',[],'KP_deformed',[],'B_deformed',[],'x_initial',[],'conn',[],'x_current',[],...
    'nnp_new',[],'u_new',[],'pn',[],'XI_orig',[],'ETA_orig',[],'elem',[],'Int_np_struc',[],'Int_Dofs',[],'Extract_op_struc',[],'Extract_op',[],...
    'K_m_struc',[],'K_m',[],'E_4',[]);


patches = repmat(patchstruc,np,1);



%defining the patches
%patch 1
patch1.pn = 1;
patch1.KP = [   0   5
    1.6     3
    3.2     1
    5       5
    5       3
    5       1
    10      5
    8.4     3
    6.8     1];
patch1.w8 = [1;1;1;
    1;1;1;
    1;1;1;];
patch1.ndm = size(patch1.KP,2);
patch1.tdm = patch1.ndm;
patch1.XI = [0 0 0  1 1 1];
patch1.ETA = [0 0 0  1 1 1];
patch1.XI_orig = [];
patch1.ETA_orig = [];
[patch1.nkn_XI, patch1.nekn_XI,  patch1.XI_elem] = element_extraction(patch1.XI);
[patch1.nkn_ETA,patch1.nekn_ETA,patch1.ETA_elem] = element_extraction(patch1.ETA);
patch1.poly_degree = 2;
patch1.p = patch1.poly_degree;
patch1.q = patch1.p;
patch1.n = patch1.nkn_XI - patch1.p - 1;
patch1.m = patch1.nkn_ETA -patch1.q - 1;
%Interface matrix
%1 number of the interface, 
%2 Patch number of the meeting patch,
%3 Master patch of the interface
%4 End basis function at the interface (0=n, 1 = m)
%5 Slave patch
patch1.Int =[1 2 0 0 1];
for i = 1:size(patch1.Int,1)
    patch1.Int(i,3) = IM(patch1.Int(i,1),4);
end
patch1.B = [];
patch1.w = [];
patch1.h_refinement_flag = [];
patch1.new_knots_XI = [];
patch1.new_knots_ETA = [];
patch1.nel = [];
patch1.nnp = [];
patch1.nen = [];
patch1.nqp = [];
patch1.ndf = [];
patch1.tdm = [];
patch1.u = [];
patch1.K = [];
patch1.fint = [];
patch1.fext = [];
patch1.fvol = [];
patch1.frea = [];
patch1.u_pre = [];
patch1.f_pre = [];
patch1.f_b = [];
patch1.IEN = [];
patch1.INC = [];
patch1.drclDofs = [];
patch1.freeDofs = [];
patch1.fpre = [];
patch1.u_d = [];
patch1.cellstruct = [];
patch1.celldata = [];
patch1.u_reshaped = [];
patch1.KP_deformed = [];
patch1.B_deformed = [];
patch1.x_initial = [];
patch1.conn = [];
patch1.x_current = [];
patch1.nnp_new = [];
patch1.u_new = [];
patch1.elem = [];
% initialize arrays for interface node points
patch1.Int_np_struc = struct('KP_n',[],'KP_c',[],'KP_f',[]);
patch1.Int_Dofs   = repmat(patch1.Int_np_struc,size(patch1.Int,1),1);
% initialize arrays for extraction operators
patch1.Extract_op_struc = struct('A1',[],'A2',[],'Csm',[],'A1_1D',[],'A2_1D',[],'Csm_1D',[]);
patch1.Extract_op   = repmat(patch1.Extract_op_struc,size(patch1.Int,1),1);
% initialize arrays for partioned stiffnes matrices
patch1.K_m_struc = struct('K_nn',[],'K_nc',[],'K_nf',[],'K_cn',[],'K_cc',[],'K_cf',[],'K_fn',[],'K_fc',[],'K_ff',[]);
patch1.K_m   = repmat(patch1.K_m_struc,size(patch1.Int,1),1);



%patch 2
patch2.pn = 2;
patch2.KP = [   3.2   1
    3.6       .5
    4       0
    5       1
    5       .5
    5       0
    6.8     1
    6.4     .5
    6       0];
    
   
patch2.w8 = [1;1;1;
    1;1;1;
    1;1;1;];
patch2.ndm = size(patch2.KP,2);
patch2.tdm = patch2.ndm;
patch2.XI = [0 0 0  1 1 1];
patch2.ETA = [0 0 0 1 1 1];
patch2.XI_orig = [];
patch2.ETA_orig = [];
[patch2.nkn_XI, patch2.nekn_XI,  patch2.XI_elem] = element_extraction(patch2.XI);
[patch2.nkn_ETA,patch2.nekn_ETA,patch2.ETA_elem] = element_extraction(patch2.ETA);
patch2.poly_degree = 2;
patch2.p = patch2.poly_degree;
patch2.q = patch2.p;
patch2.n = patch2.nkn_XI - patch2.p - 1;
patch2.m = patch2.nkn_ETA -patch2.q - 1;
%Interface matrix
%1 number of the interface, 
%2 Patch number of the meeting patch,
%3 Master patch of the interface
%4 End basis function at the interface (0=n, 1 = m)
%5 Slave patch
patch2.Int =[1 1 0 0 1];
for i = 1:size(patch2.Int,1)
    patch2.Int(i,3) = IM(patch2.Int(i,1),4);
end
patch2.B = [];
patch2.w = [];
patch2.h_refinement_flag = [];
patch2.new_knots_XI = [];
patch2.new_knots_ETA = [];
patch2.nel = [];
patch2.nnp = [];
patch2.nen = [];
patch2.nqp = [];
patch2.ndf = [];
patch2.tdm = [];
patch2.u = [];
patch2.K = [];
patch2.fint = [];
patch2.fext = [];
patch2.fvol = [];
patch2.frea = [];
patch2.u_pre = [];
patch2.f_pre = [];
patch2.f_b = [];
patch2.IEN = [];
patch2.INC = [];
patch2.drclDofs = [];
patch2.freeDofs = [];
patch2.fpre = [];
patch2.u_d = [];
patch2.cellstruct = [];
patch2.celldata = [];
patch2.u_reshaped = [];
patch2.KP_deformed = [];
patch2.B_deformed = [];
patch2.x_initial = [];
patch2.conn = [];
patch2.x_current = [];
patch2.nnp_new = [];
patch2.u_new = [];
patch2.elem = [];
% initialize arrays for interface node points
patch2.Int_np_struc = struct('KP_n',[],'KP_c',[],'KP_f',[]);
patch2.Int_Dofs   = repmat(patch2.Int_np_struc,size(patch2.Int,1),1);
% initialize arrays for extraction operators
patch2.Extract_op_struc = struct('A1',[],'A2',[],'Csm',[],'A1_1D',[],'A2_1D',[],'Csm_1D',[]);
patch2.Extract_op   = repmat(patch2.Extract_op_struc,size(patch2.Int,1),1);
% initialize arrays for partioned stiffnes matrices
patch2.K_m_struc = struct('K_nn',[],'K_nc',[],'K_nf',[],'K_cn',[],'K_cc',[],'K_cf',[],'K_fn',[],'K_fc',[],'K_ff',[]);
patch2.K_m   = repmat(patch2.K_m_struc,size(patch2.Int,1),1);

%write the patches into the patch structure
patches(1,:) = patch1;
patches(2,:) = patch2;





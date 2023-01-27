% Input.m for 2D truss structure
close all;
%clear all;
clc;
 
 %%% geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mesh = 4;

if mesh==1
     % coordinates of control points Bi,j
    %   x           y           
    %L/H = 3 fixed
      KP= [0        0           
         5          0
         10          0 
         15          0           
         20          0
         0          .5
         5          .5          
         10          .5
         15          .5
         20          .5
         0          1
         5          1          
         10          1
         15          1
         20          1
         ];

    % extract number of dsimensions from the geometry
    % defined by the control points KP
    ndm = size(KP,2);
    tdm = ndm;

    % weights 
    w8= [1;1;1;                              
         1;1;1;
         1;1;1;
         1;1;1;
         1;1;1];                            

    % % knot vectors
    XI = [0 0 0 0 1 1 1 1];
    ETA = [0 0 0 1 1 1];    

    % define the polynomial degree of the basis functions
    poly_degree = 2;

elseif mesh==2
    % coordinates of control points Bi,j
    %   x           y           
    %L/H = 3 fixed
      KP= [0        0           
         20/3       0
         40/3       0
         20         0
         0          2/3           
         20/3       2/3
         40/3       2/3
         20         2/3
         0          4/3
         20/3       4/3
         40/3       4/3
         20         4/3
         0          2          
         20/3       2
         40/3       2
         20         2
         ];

    % extract number of simensions from the geometry
    % defined by the control points KP
    ndm = size(KP,2);
    tdm = ndm;

    % weights 
    w8= [1;1;1;1;                              
         1;1;1;1;
         1;1;1;1;
         1;1;1;1];                            

    % % knot vectors
    XI = [0 0 0 0 1 1 1 1];
    ETA = [0 0 0 0 1 1 1 1];    

    % define the polynomial degree of the basis functions
    poly_degree = 3;
    
 elseif mesh==3
    % coordinates of control points Bi,j
    %   x           y           
    %L/H = 3 fixed
      KP= [0        0           
         20/4        0
         40/4       0
         60/4       0
         20          0
         0          2/4           
         20/4        2/4
         40/4       2/4
         60/4       2/4
         20          2/4
         0          1          
         20/4        1
         40/4       1
         60/4       1
         20          1
         0          6/4         
         20/4        6/4
         40/4       6/4
         60/4       6/4
         20          6/4
         0          2         
         20/4        2
         40/4       2
         60/4       2
         20          2
         ];

    % extract number of simensions from the geometry
    % defined by the control points KP
    ndm = size(KP,2);
    tdm = ndm;

    % weights 
    w8= [1;1;1;1;1;                             
         1;1;1;1;1;
         1;1;1;1;1;
         1;1;1;1;1;
         1;1;1;1;1];                            

    % % knot vectors
    XI = [0 0 0 0 0 1 1 1 1 1];
    ETA = [0 0 0 0 0 1 1 1 1 1];    

    % define the polynomial degree of the basis functions
    poly_degree = 4;
    
elseif mesh ==4
    np = 2; %number of patches
    
    %interface matrix 1: number of the index 2%3: number of the meeting
    %patches 4: master patch of the interface 5 Interface inside the
    %interface master patch matrix 6:first coloum inside the global stiffnes
    %matrix 7: last coloum inside the global stiffnes matrix
    IM = [1,1,2,1,1,0,0];
          %2,2,3,3,1,0,0];
    %[IM,index] = sortrows(IM,4:5);
    %Patch matrix 1: patch number 2:begin inside the global stiffnes matrix
    %3: end inside the global stiffnes matrix
    PM =[1,0,0;
        2,0,0];
        %3,0,0];
    CM = [1,0,0];
        %2,0,0];
    
    
    %initialize array for patch information
    patchstruc = struct('KP',[],'w8',[],'Int',[],'ndm',[],'XI',[],'ETA',[],'nkn_XI',[],'nekn_XI',[],'XI_elem',[],'nkn_ETA',[],'nekn_ETA',[],'ETA_elem',[],'poly_degree',[],'p',[],'q',[],'n',[],'m',[],...
                       'B',[],'w',[],'h_refinement_flag',[],'new_knots_XI',[],'new_knots_ETA',[],'nel',[],'nnp',[],'nen',[],'nqp',[],'ndf',[],'tdm',[],'u',[],'u_n',[],'u_f',[],'K',[],...
                       'fint',[],'fext',[],'fvol',[],'frea',[],'KP_f',[],'KP_n',[],'KP_Interface',[],'u_pre',[],'f_pre',[],'f_b',[],'IEN',[],'INC',[],'drclDofs',[],'freeDofs',[],'fpre',[],...
                   'u_d',[],'cellstruct',[],'celldata',[],'K_nn',[],'K_nf',[],'K_fn',[],'K_ff',[],'u_reshaped',[],'KP_deformed',[],'B_deformed',[],'x_initial',[],'conn',[],'x_current',[],...
                   'nnp_new',[],'u_new',[],'Ph',[],'IM',[],'pn',[],'KP_c',[],'XI_orig',[],'ETA_orig',[],'T',[],'elem',[],'Int_np_struc',[],'Int_Dofs',[],'Extract_op_struc',[],'Extract_op',[],...
                   'K_m_struc',[],'K_m',[]);
                   
    patches = repmat(patchstruc,np,1);
    
    
    
    %defining the patches
    %patch 1
    patch1.pn = 1;
    patch1.KP = [   0   0
                    5  0
                    10  0
                    0   1
                    5 1
                    10  1
                    0   2
                    5  2
                    10  2];           
    patch1.w8 = [1;1;1;
                 1;1;1;
                 1;1;1;];
    patch1.ndm = size(patch1.KP,2);
    patch1.tdm = patch1.ndm;
    patch1.XI = [0 0 0  1 1 1];
    patch1.ETA = [0 0 0  1 1 1];
    patch1.XI_orig = [];
    patch1.ETA_orig = [];
    patch1.T = [];
    [patch1.nkn_XI, patch1.nekn_XI,  patch1.XI_elem] = element_extraction(patch1.XI);
    [patch1.nkn_ETA,patch1.nekn_ETA,patch1.ETA_elem] = element_extraction(patch1.ETA);
    patch1.poly_degree = 2;
    patch1.p = patch1.poly_degree;
    patch1.q = patch1.p;
    patch1.n = patch1.nkn_XI - patch1.p - 1;
    patch1.m = patch1.nkn_ETA -patch1.q - 1;
    %Interface matrices
    %1 number of the interface, %2 x-coordinate interface start, %3 y-coordinate interface start,
    %4 x-coordinate interface end, %5 y-coordinate interface end, 
    %6 common coordinate, %7 Index common coordinate(x=1, y=2), 
    %8 Not index of the common coordinate(x=1, y=2), %9form of the
    %interface (0=abzisse m = inf, 1= straight line) 
    %10 Patch number of the meeting patch, 
    %11 Master patch of the interface 
    %12 End basis function at the interface (0=n, 1 = m)
    %13 Position of the Interface in terms of the basis function number...
    %(1-n or 1-m)
    %14 Slave patch
    patch1.Int =[1 10 0 10 2 10 1 1 0 2 0 0 patch1.n 2];
    for i = 1:size(patch1.Int,1)
        patch1.Int(i,11) = IM(patch1.Int(i,1),4);
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
    patch1.u_n = [];
    patch1.u_f = [];
    patch1.K = [];
    patch1.fint = [];
    patch1.fext = [];
    patch1.fvol = [];
    patch1.frea = [];
    patch1.KP_f = [];
    patch1.KP_c = [];
    patch1.KP_n = [];
    patch1.KP_Interface = [];
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
    patch1.K_nn = [];
    patch1.K_nf = [];
    patch1.K_fn = [];
    patch1.K_ff = [];
    patch1.u_reshaped = [];
    patch1.KP_deformed = [];
    patch1.B_deformed = [];
    patch1.x_initial = [];
    patch1.conn = [];
    patch1.x_current = [];
    patch1.nnp_new = [];
    patch1.u_new = [];
    patch1.elem = [];
    patch1.Ph = [];
    patch1.IM =[0, 0;
                1, 2;
                0, 0;
                0, 0];
    % initialize arrays for interface node points
    patch1.Int_np_struc = struct('KP_n',[],'KP_c',[],'KP_f',[]);
    patch1.Int_Dofs   = repmat(patch1.Int_np_struc,size(patch1.Int,1),1);
    % initialize arrays for extraction operators
    patch1.Extract_op_struc = struct('A1',[],'A2',[],'Csm',[]);
    patch1.Extract_op   = repmat(patch1.Extract_op_struc,size(patch1.Int,1),1);
    % initialize arrays for partioned stiffnes matrices
    patch1.K_m_struc = struct('K_nn',[],'K_nc',[],'K_nf',[],'K_cn',[],'K_cc',[],'K_cf',[],'K_fn',[],'K_fc',[],'K_ff',[]);
    patch1.K_m   = repmat(patch1.K_m_struc,size(patch1.Int,1),1);
       
                
    
    %patch 2
    patch2.pn = 2;
    patch2.KP = [   10   0
                    15  0
                    20  0
                    10   1
                    15  1
                    20  1
                    10   2
                    15  2
                    20  2];             
    patch2.w8 = [1;1;1;
                 1;1;1;
                 1;1;1;];
    patch2.ndm = size(patch2.KP,2);
    patch2.tdm = patch2.ndm;
    patch2.XI = [0 0 0  1 1 1];
    patch2.ETA = [0 0 0 1 1 1];
    patch2.XI_orig = [];
    patch2.ETA_orig = [];
    patch2.T = [];
    [patch2.nkn_XI, patch2.nekn_XI,  patch2.XI_elem] = element_extraction(patch2.XI);
    [patch2.nkn_ETA,patch2.nekn_ETA,patch2.ETA_elem] = element_extraction(patch2.ETA);
    patch2.poly_degree = 2;
    patch2.p = patch2.poly_degree;
    patch2.q = patch2.p;
    patch2.n = patch2.nkn_XI - patch2.p - 1;
    patch2.m = patch2.nkn_ETA -patch2.q - 1;
    patch2.Int =[1 10 0 10 2 10 1 2 0 1 0 0 1 2];
        %2 20 0 20 2 20 1 2 0 3 0 0 patch1.n 2];
    for i = 1:size(patch2.Int,1)
        patch2.Int(i,11) = IM(patch2.Int(i,1),4);
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
    patch2.u_n = [];
    patch2.u_f = [];
    patch2.K = [];
    patch2.fint = [];
    patch2.fext = [];
    patch2.fvol = [];
    patch2.frea = [];
    patch2.KP_f = [];
    patch2.KP_c = [];
    patch2.KP_n = [];
    patch2.KP_Interface = [];
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
    patch2.K_nn = [];
    patch2.K_nf = [];
    patch2.K_fn = [];
    patch2.K_ff = [];
    patch2.u_reshaped = [];
    patch2.KP_deformed = [];
    patch2.B_deformed = [];
    patch2.x_initial = [];
    patch2.conn = [];
    patch2.x_current = [];
    patch2.nnp_new = [];
    patch2.u_new = [];
    patch2.elem = [];
    patch2.Ph = [];
    patch2.IM =[0, 0;
        0, 0;
        0, 0;
        1, 1];
    % initialize arrays for interface node points
    patch2.Int_np_struc = struct('KP_n',[],'KP_c',[],'KP_f',[]);
    patch2.Int_Dofs   = repmat(patch2.Int_np_struc,size(patch2.Int,1),1);
    % initialize arrays for extraction operators
    patch2.Extract_op_struc = struct('A1',[],'A2',[]);
    patch2.Extract_op   = repmat(patch2.Extract_op_struc,size(patch2.Int,1),1);
    % initialize arrays for partioned stiffnes matrices
    patch2.K_m_struc = struct('K_nn',[],'K_nc',[],'K_nf',[],'K_cn',[],'K_cc',[],'K_cf',[],'K_fn',[],'K_fc',[],'K_ff',[]);
    patch2.K_m   = repmat(patch2.K_m_struc,size(patch2.Int,1),1);
       
    
%     %patch 3
%     patch3.pn = 3;
%     patch3.KP = [   20   0
%                     25  0
%                     30  0
%                     20   1
%                     25 1
%                     30  1
%                     20   2
%                     25  2
%                     30  2];           
%     patch3.w8 = [1;1;1;
%                  1;1;1;
%                  1;1;1;];
%     patch3.ndm = size(patch3.KP,2);
%     patch3.tdm = patch3.ndm;
%     patch3.XI = [0 0 0  1 1 1];
%     patch3.ETA = [0 0 0  1 1 1];
%     patch3.XI_orig = [];
%     patch3.ETA_orig = [];
%     patch3.T = [];
%     [patch3.nkn_XI, patch3.nekn_XI,  patch3.XI_elem] = element_extraction(patch3.XI);
%     [patch3.nkn_ETA,patch3.nekn_ETA,patch3.ETA_elem] = element_extraction(patch3.ETA);
%     patch3.poly_degree = 2;
%     patch3.p = patch3.poly_degree;
%     patch3.q = patch3.p;
%     patch3.n = patch3.nkn_XI - patch3.p - 1;
%     patch3.m = patch3.nkn_ETA -patch3.q - 1;
%     %Interface matrices
%     %1 number of the interface, %2 x-coordinate interface start, %3 y-coordinate interface start,
%     %4 x-coordinate interface end, %5 y-coordinate interface end, 
%     %6 common coordinate, %7 Index common coordinate(x=1, y=2), 
%     %8 Not index of the common coordinate(x=1, y=2), %9form of the
%     %interface (0=abzisse m = inf, 1= straight line) 
%     %10 Patch number of the meeting patch, 
%     %11 Master patch of the interface 
%     %12 End basis function at the interface (0=n, 1 = m)
%     %13 Position of the Interface in terms of the basis function number...
%     %(1-n or 1-m)
%     %14 Slave patch
%     patch3.Int =[2 20 0 20 2 20 1 2 0 2 0 0 1 3];
%     for i = 1:size(patch3.Int,1)
%         patch3.Int(i,11) = IM(patch3.Int(i,1),4);
%     end
%     patch3.B = [];
%     patch3.w = [];
%     patch3.h_refinement_flag = [];
%     patch3.new_knots_XI = [];
%     patch3.new_knots_ETA = [];
%     patch3.nel = [];
%     patch3.nnp = [];
%     patch3.nen = [];
%     patch3.nqp = [];
%     patch3.ndf = [];
%     patch3.tdm = [];
%     patch3.u = [];
%     patch3.u_n = [];
%     patch3.u_f = [];
%     patch3.K = [];
%     patch3.fint = [];
%     patch3.fext = [];
%     patch3.fvol = [];
%     patch3.frea = [];
%     patch3.KP_f = [];
%     patch3.KP_c = [];
%     patch3.KP_n = [];
%     patch3.KP_Interface = [];
%     patch3.u_pre = [];
%     patch3.f_pre = [];
%     patch3.f_b = [];
%     patch3.IEN = [];
%     patch3.INC = [];
%     patch3.drclDofs = [];
%     patch3.freeDofs = [];
%     patch3.fpre = [];
%     patch3.u_d = [];
%     patch3.cellstruct = [];
%     patch3.celldata = [];
%     patch3.K_nn = [];
%     patch3.K_nf = [];
%     patch3.K_fn = [];
%     patch3.K_ff = [];
%     patch3.u_reshaped = [];
%     patch3.KP_deformed = [];
%     patch3.B_deformed = [];
%     patch3.x_initial = [];
%     patch3.conn = [];
%     patch3.x_current = [];
%     patch3.nnp_new = [];
%     patch3.u_new = [];
%     patch3.elem = [];
%     patch3.Ph = [];
%     patch3.IM =[0, 0;
%         1, 2;
%         0, 0;
%         0, 0];
%     % initialize arrays for interface node points
%     patch3.Int_np_struc = struct('KP_n',[],'KP_c',[],'KP_f',[]);
%     patch3.Int_Dofs   = repmat(patch3.Int_np_struc,size(patch3.Int,1),1);
%     % initialize arrays for extraction operators
%     patch3.Extract_op_struc = struct('A1',[],'A2',[]);
%     patch3.Extract_op   = repmat(patch3.Extract_op_struc,size(patch3.Int,1),1);
%     % initialize arrays for partioned stiffnes matrices
%     patch3.K_m_struc = struct('K_nn',[],'K_nc',[],'K_nf',[],'K_cn',[],'K_cc',[],'K_cf',[],'K_fn',[],'K_fc',[],'K_ff',[]);
%     patch3.K_m   = repmat(patch3.K_m_struc,size(patch3.Int,1),1);
       

    patches(1,:) = patch1;
    patches(2,:) = patch2;
%      patches(3,:) = patch3;
    
    
elseif mesh ==5
    np = 2; %number of patches
    
    %initialize array for patch information
    patchstruc = struct('KP',[],'w8',[],'ndm',[],'XI',[],'ETA',[],'nkn_XI',[],'nekn_XI',[],'XI_elem',[],'nkn_ETA',[],'nekn_ETA',[],'ETA_elem',[],'poly_degree',[],'p',[],'q',[],'n',[],'m',[],...
                       'B',[],'w',[],'h_refinement_flag',[],'new_knots_XI',[],'new_knots_ETA',[],'nel',[],'nnp',[],'nen',[],'nqp',[],'ndf',[],'tdm',[],'u',[],'u_n',[],'u_f',[],'K',[],...
                       'fint',[],'fext',[],'fvol',[],'frea',[],'KP_f',[],'KP_n',[],'KP_Interface',[],'u_pre',[],'f_pre',[],'f_b',[],'IEN',[],'INC',[],'drclDofs',[],'freeDofs',[],'fpre',[],...
                   'u_d',[],'cellstruct',[],'celldata',[],'K_nn',[],'K_nf',[],'K_fn',[],'K_ff',[],'u_reshaped',[],'KP_deformed',[],'B_deformed',[],'x_initial',[],'conn',[],'x_current',[],'nnp_new',[],'u_new',[]);
                   
    patches = repmat(patchstruc,np,1);
    
    
    
    %defining the patches
    %patch 1
    patch1.KP = [10  0
                11   0
                0  1
                1  1  ];
    patch1.w8 = [1;1;                              
                1;1;];
    patch1.ndm = size(patch1.KP,2);
    patch1.tdm = patch1.ndm;
    patch1.XI = [0 0 1 1];
    patch1.ETA = [0 0 1 1];
    [patch1.nkn_XI, patch1.nekn_XI,  patch1.XI_elem] = element_extraction(patch1.XI);
    [patch1.nkn_ETA,patch1.nekn_ETA,patch1.ETA_elem] = element_extraction(patch1.ETA);
    patch1.poly_degree = 1;
    patch1.p = patch1.poly_degree;
    patch1.q = patch1.p;
    patch1.n = patch1.nkn_XI - patch1.p - 1;
    patch1.m = patch1.nkn_ETA -patch1.q - 1;
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
    patch1.u_n = [];
    patch1.u_f = [];
    patch1.K = [];
    patch1.fint = [];
    patch1.fext = [];
    patch1.fvol = [];
    patch1.frea = [];
    patch1.KP_f = [];
    patch1.KP_n = [];
    patch1.KP_Interface = [];
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
    patch1.K_nn = [];
    patch1.K_nf = [];
    patch1.K_fn = [];
    patch1.K_ff = [];
    patch1.u_reshaped = [];
    patch1.KP_deformed = [];
    patch1.B_deformed = [];
    patch1.x_initial = [];
    patch1.conn = [];
    patch1.x_current = [];
    patch1.nnp_new = [];
    patch1.u_new = [];
    
    %patch 2
    patch2.KP = [1  0
                 2  0
                 1  1
                 2  1];
    patch2.w8 = [1;1;                              
                1;1];
    patch2.ndm = size(patch2.KP,2);
    patch2.tdm = patch2.ndm;
    patch2.XI = [0 0 1 1];
    patch2.ETA = [0 0 1 1];
    [patch2.nkn_XI, patch2.nekn_XI,  patch2.XI_elem] = element_extraction(patch2.XI);
    [patch2.nkn_ETA,patch2.nekn_ETA,patch2.ETA_elem] = element_extraction(patch2.ETA);
    patch2.poly_degree = 1;
    patch2.p = patch2.poly_degree;
    patch2.q = patch2.p;
    patch2.n = patch2.nkn_XI - patch2.p - 1;
    patch2.m = patch2.nkn_ETA -patch2.q - 1;
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
    patch2.u_n = [];
    patch2.u_f = [];
    patch2.K = [];
    patch2.fint = [];
    patch2.fext = [];
    patch2.fvol = [];
    patch2.frea = [];
    patch2.KP_f = [];
    patch2.KP_n = [];
    patch2.KP_Interface = [];
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
    patch2.K_nn = [];
    patch2.K_nf = [];
    patch2.K_fn = [];
    patch2.K_ff = [];
    patch2.u_reshaped = [];
    patch2.KP_deformed = [];
    patch2.B_deformed = [];
    patch2.x_initial = [];
    patch2.conn = [];
    patch2.x_current = [];
    patch2.nnp_new = [];
    patch2.u_new = [];

    patches(1,:) = patch1;
    patches(2,:) = patch2;
    

end



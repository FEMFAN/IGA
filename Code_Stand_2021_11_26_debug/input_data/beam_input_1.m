% Input.m for 2D truss structure
close all;
clear all;
clc;
 
 %%% geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mesh = 5;

if mesh==1
     % coordinates of control points Bi,j
    %   x           y           
    %L/H = 3 fixed
      KP= [0        0           
         5          0
         10          0 
         15          0           
         20          0
         0          1
         5          1          
         10          1
         15          1
         20          1
         0           2
         5          2
         10         2
         15         2
         20         2
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
    
    %initialize array for patch information
    patchstruc = struct('KP',[],'w8',[],'ndm',[],'XI',[],'ETA',[],'nkn_XI',[],'nekn_XI',[],'XI_elem',[],'nkn_ETA',[],'nekn_ETA',[],'ETA_elem',[],'poly_degree',[],'p',[],'q',[],'n',[],'m',[],...
                       'B',[],'w',[],'h_refinement_flag',[],'new_knots_XI',[],'new_knots_ETA',[],'nel',[],'nnp',[],'nen',[],'nqp',[],'ndf',[],'tdm',[],'u',[],'u_n',[],'u_f',[],'K',[],...
                       'fint',[],'fext',[],'fvol',[],'frea',[],'KP_f',[],'KP_n',[],'KP_Interface',[],'u_pre',[],'f_pre',[],'f_b',[],'IEN',[],'INC',[],'drclDofs',[],'freeDofs',[],'fpre',[],...
                   'u_d',[],'cellstruct',[],'celldata',[],'K_nn',[],'K_nf',[],'K_fn',[],'K_ff',[],'u_reshaped',[],'KP_deformed',[],'B_deformed',[],'x_initial',[],'conn',[],'x_current',[],'nnp_new',[],'u_new',[]);
                   
    patches = repmat(patchstruc,np,1);
    
    
    
    %defining the patches
    %patch 1
    patch1.KP = [0  0
                5   0
                10  0
                0   1
                5   1
                10  1
                0   2
                5   2
                10  2
                ];
    patch1.w8 = [1;1;1;                              
                1;1;1;
                1;1;1
                ];
    patch1.ndm = size(patch1.KP,2);
    patch1.tdm = patch1.ndm;
    patch1.XI = [0 0 0 1 1 1];
    patch1.ETA = [0 0 0 1 1 1];
    [patch1.nkn_XI, patch1.nekn_XI,  patch1.XI_elem] = element_extraction(patch1.XI);
    [patch1.nkn_ETA,patch1.nekn_ETA,patch1.ETA_elem] = element_extraction(patch1.ETA);
    patch1.poly_degree = 2;
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
    patch2.KP = [10  0
                15   0
                20  0
                10   1
                15   1
                20  1
                10   2
                15   2
                20  2];
    patch2.w8 = [1;1;1;                              
                1;1;1;
                1;1;1];
    patch2.ndm = size(patch2.KP,2);
    patch2.tdm = patch2.ndm;
    patch2.XI = [0 0 0 1 1 1];
    patch2.ETA = [0 0 0 1 1 1];
    [patch2.nkn_XI, patch2.nekn_XI,  patch2.XI_elem] = element_extraction(patch2.XI);
    [patch2.nkn_ETA,patch2.nekn_ETA,patch2.ETA_elem] = element_extraction(patch2.ETA);
    patch2.poly_degree = 2;
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
    patch1.KP = [0  0
                1   0
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



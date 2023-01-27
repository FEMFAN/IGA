function [xi_tilde,w8_tilde] = gaussian_quadrature(nqp,ndm)
% ---------------------------------------------------------------------
% Subroutine gauss2d.m
% calculates gauss point coordiantes and weights
% for numerical integration
% Source: Vinh Phu Nguyen
%
% Author:           Carina Witt
% Date  :           27.11.2018
%
% Input:    nqp             - number of integration points
%           ndm             - number of dimensions
%
% Output:   xi_tilde        - gauss point coordinates
%           w8_tilde        - integration weights
%-----------------------------------------------------------------------

% check if quadrature point value is valid
if nqp<4
    error('Gaussian quadrature for nqp not implemented!');  
end

% calculate quadrature order
quadorder = nqp^(1/ndm);

% initialize gauss points and weights for the quadrature order 
qpoint  = zeros(quadorder,1); 
qweight = zeros(quadorder,1);
            
% assign quadrature point values and their weights
% depending on the quadrature order
if quadorder==2
    qpoint(1) = 0.577350269189626;
    qpoint(2) =-0.577350269189626;

    qweight(1) = 1.000000000000000;
    qweight(2) = 1.000000000000000;
            
elseif quadorder==3
    qpoint(1) = 0.774596669241483;
    qpoint(2) =-0.774596669241483;
    qpoint(3) = 0.000000000000000;

    qweight(1) = 0.555555555555556;
    qweight(2) = 0.555555555555556;
    qweight(3) = 0.888888888888889;
            
elseif quadorder==4
    qpoint(1) = 0.861134311594053;
    qpoint(2) =-0.861134311594053;
    qpoint(3) = 0.339981043584856;
    qpoint(4) =-0.339981043584856;
            
    qweight(1) = 0.347854845137454;
    qweight(2) = 0.347854845137454;
    qweight(3) = 0.652145154862546;
    qweight(4) = 0.652145154862546;

elseif quadorder==5
    qpoint(1) = 0.906179845938664;
    qpoint(2) =-0.906179845938664;
    qpoint(3) = 0.538469310105683;
    qpoint(4) =-0.538469310105683;
    qpoint(5) = 0.000000000000000;

    qweight(1) = 0.236926885056189;
    qweight(2) = 0.236926885056189;
    qweight(3) = 0.478628670499366;
    qweight(4) = 0.478628670499366;
    qweight(5) = 0.568888888888889;

elseif quadorder==6
    qpoint(1) = 0.932469514203152;
    qpoint(2) =-0.932469514203152;
    qpoint(3) = 0.661209386466265;
    qpoint(4) =-0.661209386466265;
    qpoint(5) = 0.238619186003152;
    qpoint(6) =-0.238619186003152;

    qweight(1) = 0.171324492379170;
    qweight(2) = 0.171324492379170;
    qweight(3) = 0.360761573048139;
    qweight(4) = 0.360761573048139;
    qweight(5) = 0.467913934572691;
    qweight(6) = 0.467913934572691;

elseif quadorder==7
    qpoint(1) =  0.949107912342759;
    qpoint(2) = -0.949107912342759;
    qpoint(3) =  0.741531185599394;
    qpoint(4) = -0.741531185599394;
    qpoint(5) =  0.405845151377397;
    qpoint(6) = -0.405845151377397;
    qpoint(7) =  0.000000000000000;

    qweight(1) = 0.129484966168870;
    qweight(2) = 0.129484966168870;
    qweight(3) = 0.279705391489277;
    qweight(4) = 0.279705391489277;
    qweight(5) = 0.381830050505119;
    qweight(6) = 0.381830050505119;
    qweight(7) = 0.417959183673469;

elseif quadorder==8
    qpoint(1) =  0.960289856497536;
    qpoint(2) = -0.960289856497536;
    qpoint(3) =  0.796666477413627;
    qpoint(4) = -0.796666477413627;
    qpoint(5) =  0.525532409916329;
    qpoint(6) = -0.525532409916329;
    qpoint(7) =  0.183434642495650;
    qpoint(8) = -0.183434642495650;

    qweight(1) = 0.101228536290376;
    qweight(2) = 0.101228536290376;
    qweight(3) = 0.222381034453374;
    qweight(4) = 0.222381034453374;
    qweight(5) = 0.313706645877887;
    qweight(6) = 0.313706645877887;
    qweight(7) = 0.362683783378362;
    qweight(8) = 0.362683783378362;
elseif quadorder==9
    qpoint(1) = -0.968160239507626;
    qpoint(2) = -0.836031107326636;
    qpoint(3) = -0.613371432700590;
    qpoint(4) = -0.324253423403809;
    qpoint(5) =  0.000000000000000;
    qpoint(6) =  0.324253423403809;
    qpoint(7) =  0.613371432700590;
    qpoint(8) =  0.836031107326636;
    qpoint(9) =  0.968160239507626;

    qweight(1) =  0.081274388361574;
    qweight(2) =  0.180648160694857;
    qweight(3) =  0.260610696402935;
    qweight(4) =  0.312347077040003;
    qweight(5) =  0.330239355001260;
    qweight(6) =  0.312347077040003;
    qweight(7) =  0.261610696402935;
    qweight(8) =  0.180648160694857;
    qweight(9) =  0.081274388361574;
elseif quadorder==10
    qpoint(1)  = -0.973906528517172;
    qpoint(2)  = -0.865063366688985;
    qpoint(3)  = -0.679409568299024;
    qpoint(4)  = -0.433395394129247;
    qpoint(5)  = -0.148874338981631;
    qpoint(6)  =  0.148874338981631;
    qpoint(7)  =  0.433395394129247;
    qpoint(8)  =  0.679409568299024;
    qpoint(9)  =  0.865063366688985;
    qpoint(10) =  0.973906528517172;

    qweight(1)  =  0.066671344308688;
    qweight(2)  =  0.149451349150581;
    qweight(3)  =  0.219086362515982;
    qweight(4)  =  0.269266719309996;
    qweight(5)  =  0.295524224714753;
    qweight(6)  =  0.295524224714753;
    qweight(7)  =  0.269266719309996;
    qweight(8)  =  0.219086362515982;
    qweight(9)  =  0.149451349150581;
    qweight(10) =  0.066671344308688;
elseif quadorder==12
    qpoint(1)  = -0.981560634246719;
    qpoint(2)  = -0.904117256370475;
    qpoint(3)  = -0.769902674194305;
    qpoint(4)  = -0.587317954286617;
    qpoint(5)  = -0.367831498998180;
    qpoint(6)  = -0.125233408511469;
    qpoint(7)  =  0.125233408511469;
    qpoint(8)  =  0.367831498998180;
    qpoint(9)  =  0.587317954286617;
    qpoint(10) =  0.769902674194305;
    qpoint(11) =  0.904117256370475;
    qpoint(12) =  0.981560634246719;

    qweight(1)  =  0.047175336386512;
    qweight(2)  =  0.106939325995318;
    qweight(3)  =  0.160078328543346;
    qweight(4)  =  0.203167426723066;
    qweight(5)  =  0.233492536538355;
    qweight(6)  =  0.249147045813403;
    qweight(7)  =  0.249147045813403;
    qweight(8)  =  0.233492536538355;
    qweight(9)  =  0.203167426723066;
    qweight(10) =  0.160078328543346;
    qweight(11) =  0.106939325995318;
    qweight(12) =  0.047175336386512;

else
     error('Gaussian quadrature for nqp not implemented!');  
end  % end if
   
% assemble matrix xi_tilde with gauss point values
% and vector w8_tilde with the corresponding weigths
n=1;
if ( ndm == 1 )
    for i = 1:quadorder
        xi_tilde(:,n) = [ xi_tilde(i) ];
        w8_tilde(1,n)  = qweight(i);
        n = n+1;
    end

elseif ( ndm == 2 )
    for i = 1:quadorder
        for j = 1:quadorder
            xi_tilde(:,n) = [ qpoint(i); qpoint(j)];
            w8_tilde(1,n)  = qweight(i)*qweight(j);
            n = n+1;
        end
    end

else % ndm == 3
    for i = 1:quadorder
        for j = 1:quadorder
            for k = 1:quadorder
                xi_tilde(:,n) = [ qpoint(i); qpoint(j); qpoint(k) ];
                w8_tilde(1,n) = qweight(i)*qweight(j)*qweight(k);
                n = n+1;
            end
        end
    end

end
%xi_tilde
%w8_tilde

end % function
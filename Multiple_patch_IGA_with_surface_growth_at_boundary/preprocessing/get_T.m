function T = get_T(XI_orig,XI_refined,p)
% ---------------------------------------------------------------------
% Subroutine get_T.m
% calculation of the extension operator as shown in literature Hughes and Cotrell <<ISOGEOMETRIC ANALYSIS TOWARD INTEGRATION OF CAD AND FEA>>
%
% Author:           Lennard Langerbein
% Date  :           15.02.2022
%                  
% Input:    XI_orig                             - original knot vector
%           XI_refined                          - extended knot vector
%           p                                   - polynomial degree
%
% Output:   T                                   - extension operator  
%---------------------------------------------------------------------- 
Xi_Bar = XI_refined;
Xi = XI_orig;
T0 = zeros(length(Xi_Bar),length(Xi)-p);
for i = 1:length(Xi_Bar)
    for j = 1:length(Xi)-1
        if Xi_Bar(i) >= Xi(j) && Xi_Bar(i) < Xi(j+1)
            T0(i,j) = 1;
        else
            T0(i,j) = 0;
        end
    end
end
%extension operator source Hughes and Cotrell <<ISOGEOMETRIC ANALYSIS TOWARD INTEGRATION OF CAD AND FEA>>
T0 = T0(1:end-(p+1),:);
for q = 0:(p-1)
    if q == 0
        for i = 1:size(T0,1)
            for j = 1:size(T0,2)-1
                if (Xi_Bar(i+q+1)-Xi(j)) < 1e-8 || (Xi(j+q+1)-Xi(j)) < 1e-8
                    A = 0;
                    if (Xi(j+q+1+1)-Xi_Bar(i+q+1)) < 1e-8 ||(Xi(j+q+1+1)-Xi(j+1)) < 1e-8
                        B = 0;
                    else
                        B = (Xi(j+q+1+1)-Xi_Bar(i+q+1))/(Xi(j+q+1+1)-Xi(j+1))* T0(i,j+1);
                    end
                elseif (Xi(j+q+1+1)-Xi_Bar(i+q+1)) < 1e-8 || (Xi(j+q+1+1)-Xi(j+1)) < 1e-8
                    B = 0;
                    if (Xi_Bar(i+q+1)-Xi(j)) < 1e-8 || (Xi(j+q+1)-Xi(j)) < 1e-8
                        A = 0;
                    else
                        A = (Xi_Bar(i+q+1)-Xi(j))/(Xi(j+q+1)-Xi(j))* T0(i,j);
                    end
                else
                    A = (Xi_Bar(i+q+1)-Xi(j))/(Xi(j+q+1)-Xi(j))* T0(i,j);
                    B = (Xi(j+q+1+1)-Xi_Bar(i+q+1))/(Xi(j+q+1+1)-Xi(j+1))* T0(i,j+1);
                end
                T(i,j,q+1) = A + B;
            end
        end
    else
        for i = 1:size(T(:,:,q),1)
            for j = 1:size(T(:,:,q),2)-1
                if (Xi_Bar(i+q+1)-Xi(j)) < 1e-8 || (Xi(j+q+1)-Xi(j)) < 1e-8
                    A = 0;
                    if (Xi(j+q+1+1)-Xi_Bar(i+q+1)) < 1e-8 || (Xi(j+q+1+1)-Xi(j+1)) < 1e-8
                        B = 0;
                    else
                        B = (Xi(j+q+1+1)-Xi_Bar(i+q+1))/(Xi(j+q+1+1)-Xi(j+1))* T(i,j+1,q);
                    end
                elseif (Xi(j+q+1+1)-Xi_Bar(i+q+1)) < 1e-8 || (Xi(j+q+1+1)-Xi(j+1)) < 1e-8
                    B = 0;
                    if (Xi_Bar(i+q+1)-Xi(j)) < 1e-8 || (Xi(j+q+1)-Xi(j)) < 1e-8
                        A = 0;
                    else
                        A = (Xi_Bar(i+q+1)-Xi(j))/(Xi(j+q+1)-Xi(j))* T(i,j,q);
                    end
                else
                    A = (Xi_Bar(i+q+1)-Xi(j))/(Xi(j+q+1)-Xi(j))* T(i,j,q);
                    B = (Xi(j+q+1+1)-Xi_Bar(i+q+1))/(Xi(j+q+1+1)-Xi(j+1))* T(i,j+1,q);
                end
                T(i,j,q+1) = A + B;
            end
        end
    end
end
if p==1
    T = T(:,1:end,p);
else
    T = T(:,1:end-1,p);
end
end
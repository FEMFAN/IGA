function [E4] = Hooke_material(ndm,k,material_flag)

% Hooke, isotropic elasticity 

% material parameters
if k ==2 & material_flag == 1
    xE  = .7e5;
    xnu = 0.35;
else
    xE  = 2.1e5;
    xnu = 0.3;
end


% calculate Lame parameters
xlm = xE * xnu / ( (1.0 + xnu)*(1.0 - 2.0*xnu) );
xmu = xE / ( 2.0*(1.0 + xnu) );

% Kronecker delta
I = eye(ndm,ndm);

% 4th-order elasticity tensor
E4 = zeros(ndm,ndm,ndm,ndm);
for i=1:ndm
    for j=1:ndm
        for k=1:ndm
            for l=1:ndm
                E4(i,j,k,l) = xlm*I(i,j)*I(k,l) + 2.0*xmu*( 0.5*( I(i,k)*I(j,l) + I(i,l)*I(j,k) ) );
            end
        end
    end
end

end % function
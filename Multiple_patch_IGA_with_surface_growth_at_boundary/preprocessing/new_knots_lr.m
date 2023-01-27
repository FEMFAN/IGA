function [patches] = new_knots_lr(patches,lflag)
% ---------------------------------------------------------------------
% Subroutine new_knots_lr.m
% calculate new knots for local refinment
%
% Author:           Lennard Langerbein
% Date  :           15.02.2022
%                  
%
% Input:    patches                             -structure with the patch information
%           lflag                               - local refinement flag
%
% Output:   patches                             - structure with refined
%                                                   knot vecors             
%---------------------------------------------------------------------- 
%save the old knot vectors and compute the new ones
Xi = patches(1).XI;
Eta = patches(1).ETA;
new_knots_Xi = [];
new_knots_Eta = [];
for i = 1:length(Xi)-1
    if abs(Xi(i) - Xi(i+1)) > 1e-8
        delta = abs((Xi(i)-Xi(i+1))/lflag);
        for l = 1:lflag-1
            new_knots_Xi(end+1) = Xi(i)+l*delta;
        end
    end
end

patches(1).new_knots_XI = new_knots_Xi;
for i = 1:length(Eta)-1
    if abs(Eta(i) - Eta(i+1)) > 1e-8
        delta = abs((Eta(i+1)-Eta(i))/lflag);
        for l = 1:lflag-1
            new_knots_Eta(end+1) =Eta(i)+l*delta;
        end
    end
end
%save the new knot vector inside the patch structure
patches(1).new_knots_ETA = new_knots_Eta;
end


function [XI_mod,nkn_XI_mod,nekn_XI_mod,XI_elem_mod,B_mod,w_mod,n_mod] = knot_removal_XI(m,ndm,XI_orig,poly_degree,B,w,remove_knot)
% ---------------------------------------------------------------------
% Subroutine knot_removal_XI.m
% Defines new knot vector and control points
% on the basis of the knots to be removed
% literature: Piegl,Tiller 5.4
%
% Attention: only correct for removal multiplicity of 1
%
% Author:           M.Sc. Carina Witt
% Date  :           23.03.2020
%
% Input:    m               - number of shape functions
%           ndm             - number of dimensions
%           XI_orig         - original knot vector
%           xi              - current coordinate value
%           nkn             - original number of knots
%           poly_degree     - ploynomial degree of the basis functions
%           B               - original control points
%           w               - original weights of the control points
%           remove_knot     - knot to be removed
%
% Output:   XI_mod          - refined knot vector
%           nkn_XI_mod      - new number of knots
%           nekn_XI_mod     - new number of elements (knot spans)
%           XI_elem_mod     - new elements (knot spans)
%           KP_mod          - control points for new knot vector
%           w8_mod          - weights
%           n_mod           - new number of basis functions 
%---------------------------------------------------------------------- 

% loop through all curves 
for m_i = 1:m


    XI = XI_orig;

    % Bring controlpoint vector in correct representation for knot removal
    % knot removal is not carried out on projected control points 
    % but on the 'hyperplane'
    KPw = B(:,:,m_i).*w(:,1,m_i);
    KPw(:,end+1) = w(:,1,m_i);

    %% Calculate new control points and weights
    % before check if knot removal is applicable without changing the curve shape

for knot_i = 1:size(remove_knot,1)
    
    nkn  = size(XI,2);
    n = nkn-poly_degree-1;

    % get knot value that shall be removed and how many times
    remove_multiplicity = remove_knot(knot_i,2);
    if remove_multiplicity >1
		fprintf('Warning: knot removal might be inaccurate!')
	end
    remove_knot_ii      = remove_knot(knot_i,1);

    % initial multiplicity of the knot that shall be removed
    initial_multiplicity = 0;
    for i=2:(nkn-1) % end knots are not being removed

      % check position of the knot to be removed
      if remove_knot_ii == XI(1,i) 
        % save position of the new knot as r 
        % counter starts at zero!
        r = i-1;
        % compute initial multiplicity
        initial_multiplicity = initial_multiplicity +1;
      end

    end
    s = initial_multiplicity;

    % variables to determine locations in control point vector
    ord = poly_degree+1;
    fout = 0.5*(2*r-s-poly_degree);

    % boundaries for modification area
    % -> control points from "first" to "last" will be modified 
    first = r-poly_degree+1; last = r-s+1;

    TOL = 1;
    for t = 0:(remove_multiplicity-1)

        % Index difference between original and temporary modified control polygon
        off = first-1; 

        % initialize modified control point values
        tempKPw = zeros(2*poly_degree+1,3); % HERE : initialize before loop?!
        tempKPw(1,:)          = KPw(off,:);
        tempKPw(last+2-off,:) = KPw(last+1,:);

        % start from the outside and then move towards interior
        i = first; j = last;
        ii = 1; jj = last+1-off;

        % flag determines whether removal is applicable
        % 0 - no removal
        % 1 - removal
        remflag = 0;

        while ( (j-i)>t )

            alpha_i = ( remove_knot_ii-XI(i) )  /( XI(i+ord+t)-XI(i) );
            alpha_j = ( remove_knot_ii-XI(j-t) )/( XI(j+ord)-XI(j-t) );

            tempKPw(ii+1,:) = ( KPw(i,:) - (1-alpha_i)*tempKPw(ii,:) ) / alpha_i;
            tempKPw(jj,:)   = ( KPw(j,:) -    alpha_j *tempKPw(jj+1,:) ) / ( 1- alpha_j );
            
            i  = i+1;   j = j-1;
            ii = ii+1; jj = jj-1;

        end % while

        % check if knot is removable
        if (j-i)<t 
            if norm( tempKPw(ii,:) - tempKPw(jj+1,:) ) <= TOL % check eq. (5.26)
                % Check weight values (these can be negative with this
                % algorithm) 
                for w_i = 1:size(KPw,1)
                    if KPw(w_i,ndm+1) <= 0
                        fprintf('Warning: Negative weights in knot removal procedure!')
                    end
                end
                remflag = 1;
            end
        else       
            alpha_i = ( remove_knot_ii-XI(i) )  /( XI(i+ord+t)-XI(i) );
            if norm( KPw(i,:) - ( alpha_i*tempKPw(ii+t+2,:)+(1-alpha_i)*tempKPw(ii,:) ) ) <= TOL % check eq. (5.25), (5.27)
                % Check weight values (these can be negative with this
                % algorithm) 
                for w_i = 1:size(KPw,1)
                    if KPw(w_i,ndm+1) <= 0
                        fprintf('Warning: Negative weights in knot removal procedure!')
                    end
                end
                remflag = 1;
            end
        end

        if remflag == 0
            fprintf('Warning: Knot cannot be removed (anymore)!')
            break; % no knots to be removed
        else
            % removal successful 
            % save new control points and weights
            i = first; j = last;
            while (j-i)>t
               KPw(i,:) = tempKPw(i+1-off,:);
               KPw(j,:) = tempKPw(j+1-off,:);
               i = i+1; j = j-1;
            end
        end

        first = first-1; last = last+1;

    end % for

    % shift knot vector and control point and weight values

    % shift knot vector
    for k=(r+1):(nkn-1)
       XI(k-t) = XI(k+1);
    end
    % overwrite XI
    XI = XI(1:(end-(t+1)));

    % shift list of control points and weights
    j = fout; i=j; 
    % check even or odd case 
    if t == 0
       j = j+1; i=j;
    end

    for k=1:t 
       if mod(k,2)==1 % odd multiplicity of removal
           i=i+t; % modification on 24.03.2020 (before: i+1)
       else % even multiplicity of removal
           j=j-2+t; % modification on 24.03.2020 (before: j-1)
       end
    end
    
    for k=(i+1):n % perform shift
       KPw(j,:) = KPw(k,:);
       j = j+1;
    end
    % overwrite KPw
    KPw = KPw(1:(end-(t+1)),:);       
 
end % for

%% Obtain projected control polygon
B_mod(:,:,m_i) = KPw(:,1:ndm);  
B_mod(:,1:ndm,m_i) = B_mod(:,1:ndm,m_i)./KPw(:,ndm+1);
w_mod(:,1,m_i) = KPw(:,ndm+1);
B_mod(:,:,m_i) = B_mod(:,1:ndm,m_i);

% modify control point values manually
% aim: avoid numerical mistakes 
% -> here: set values near to zero (error < 0.05%) to exact zero
% attention: factor might not be suitable for different mesh sizes!!
% otherwise circular shape is not depicted accurately
factor = 0.0005;
for test_i = 1:size(B_mod,1)
    for ndm_i = 1:ndm
       if ( B_mod(test_i,ndm_i,m_i) > 0.0 && B_mod(test_i,ndm_i,m_i) < (factor*B_mod(1,1,m_i)) ) || ...
               ( B_mod(test_i,ndm_i,m_i) > -(factor*B_mod(1,1,m_i)) && B_mod(test_i,ndm_i,m_i) < 0.0 )
           B_mod(test_i,ndm_i,m_i) = 0.0;
       end
    end
end
KPw(:,1:ndm) = B_mod(:,:,m_i).*w_mod(:,1,m_i);
KPw(:,ndm+1) = w_mod(:,1,m_i);
            
end % loop over m 

% extract new number of knots and knot spans from knot vector
% via subroutine element_extraction.m
XI_mod = XI;
[nkn_XI_mod,nekn_XI_mod,XI_elem_mod] = element_extraction(XI_mod);
n_mod = nkn_XI_mod-poly_degree-1;

end % function

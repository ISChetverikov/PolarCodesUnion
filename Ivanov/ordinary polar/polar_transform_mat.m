function x = polar_transform_mat(u)
    % Recurse down to length 1
    if (size(u,2)==1)
        x = u;
    else
        % Compute odd/even outputs of (I_{N/2} \otimes G_2) transform
        u1u2 = mod(u(:,1:2:end)+u(:,2:2:end),2);
        u2 = u(:,2:2:end);
        % R_N maps odd/even indices (i.e., u1u2/u2) to first/second half
        x = [polar_transform_mat(u1u2) polar_transform_mat(u2)];
    end
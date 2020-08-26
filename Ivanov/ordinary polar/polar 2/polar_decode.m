function [u,x] = polar_decode(y,f)
    % y = bit APP from channel in output order % f = input a priori probs in input order
    % x = output hard decision in output order % u = input hard decisions in input order
    % Recurse down to length 1
    N = length(y);
    if (N==1)
        if (f==1/2)
            % If info bit, make hard decision based on observation
            x = (1-sign(1-2*y))/2; u = x;
        else
            % If frozen, use frozen bit
            x = f; u = x;
        end
    else
        % Compute soft mapping back one stage
        u1est = cnop(y(1:2:end),y(2:2:end));
        % R_N^T maps u1est to top polar code
        [uhat1,u1hardprev] = polar_decode(u1est,f(1:(N/2)));
        % Using u1est and x1hard, we can estimate u2
        u2est = vnop(cnop(u1hardprev,y(1:2:end)),y(2:2:end));
        % R_N^T maps u2est to bottom polar code
        [uhat2,u2hardprev] = polar_decode(u2est,f((N/2+1):end));
        % Tunnel u decisions back up. Compute and interleave x1,x2 hard decisions
        u = [uhat1 uhat2];
        x = reshape([cnop(u1hardprev,u2hardprev); u2hardprev],1,[]);
    end
    return
    % Check-node operation in P1 domain
function z=cnop(w1,w2)
    z = w1.*(1-w2) + w2.*(1-w1);
    return
    % Bit-node operation in P1 domain
function z=vnop(w1,w2)
    z = w1.*w2 ./ (w1.*w2 + (1-w1).*(1-w2));
    return


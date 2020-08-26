function [u,x,idx,listprob] = polar_decode_pcscl(y,f,p,L,listprob,u)
    % y = bit APP from channel in output order (#row = index in list)
    % f = input a priori probs in input order
    % L = maximum list size
    % listprob = probabilities of each path in the list
    % x = output hard decision in output order
    % u = input hard decisions in input order
    % idx = indices of the list items used
    % Recurse down to length 1
    if nargin < 6, u = []; end
    if nargin < 5, listprob = 1; end
    if nargin < 4, L = 1; end
    N = size(y,2);
    if (N==1)
        L0 = size(y,1);
        if f == 'i'
            % Make list of all possible decisions (double the list size)
            [listprob,j] = sort([listprob.*y; listprob.*(1-y)],1,'descend');
            % Trim list if it is too large
            if length(j) > L
                listprob = listprob(1:L);
                j = j(1:L);
            end
            % the first half of the list was due to x=1
            x = j <= L0;
            idx = mod(j-1,L0)+1;
            u = [u(idx,:) x];
            % determine which path each decision belongs to
        elseif f == 'p'
            % If pc-frozen, calculate bit value
            i = mod(size(u,2),p)+1;
            x = mod(sum(u(:,i:p:end),2),2);
            u = [u x];
            idx = 1:L0;
            listprob = listprob .* (x.*y + (1-x).*(1-y)); % Or better avoid modifing it?
        else
            % If frozen, use frozen bit
            x = zeros(size(y));
            u = [u x];
            idx = 1:L0;
            listprob = listprob .* (1-y); % Or better avoid modifing it?
        end
        listprob = listprob / max(listprob);
    else
        % Compute soft mapping back one stage
        u1est = cnop(y(:,1:2:end),y(:,2:2:end));
        % R_N^T maps u1est to top polar code
        [uhat1,u1hardprev,idx1,prob1] = polar_decode_pcscl(u1est,f(1:(N/2)),p,L,listprob,u);
        % Using u1est and x1hard, we can estimate u2
        u2est = vnop(cnop(u1hardprev,y(idx1,1:2:end)),y(idx1,2:2:end));
        % R_N^T maps u2est to bottom polar code
        [uhat2,u2hardprev,idx2,prob2] = polar_decode_pcscl(u2est,f((N/2+1):end),p,L,prob1,uhat1);
        % Tunnel u decisions back up. Compute and interleave x1,x2 hard decisions
        u = uhat2;
        x = zeros(length(idx2), size(u2hardprev,2)*2);
        x(:,1:2:end) = cnop(u1hardprev(idx2,:),u2hardprev);
        x(:,2:2:end) = u2hardprev;
        listprob = prob2;
        idx = idx1(idx2);
    end
    return
    % Check-node operation in P1 domain
function z=cnop(w1,w2)
    z = w1.*(1-w2) + w2.*(1-w1);
    return
    % Bit-node operation in P1 domain
function z=vnop(w1,w2)
    z = w1.*w2 ./ (w1.*w2 + (1-w1).*(1-w2));
    z(isnan(z)) = .5;
    return


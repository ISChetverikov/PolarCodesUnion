function [u,x,idx,listprob] = polar_decode_list(y,f,L,listprob,u)
    % y = bit APP from channel in output order (#row = index in list)
    % f = input a priori probs in input order
    % L = maximum list size
    % listprob = probabilities of each path in the list
    % x = output hard decision in output order
    % u = input hard decisions in input order
    % idx = indices of the list items used
    % Recurse down to length 1
    if nargin < 5, u = []; end
    if nargin < 4, listprob = 0; end
    if nargin < 3, L = 1; end
    if ~iscell(f), f = num2cell(f); end
    N = size(y,2);
    if (N==1)
        Lp = llr2logp(y);
        L0 = size(y,1);
        if (f{1}==1/2)
            % Make list of all possible decisions (double the list size)
            [listprob,j] = sort([listprob + Lp(:,2); listprob + Lp(:,1)],1,'descend');
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
        else
            % If frozen, use frozen bit
            %x = f; u = x;
            x = repmat(f{1},size(y));
            %x = f{1} * u;
            %if size(x,2) > 1
            %    x = sum(x,2); % if f{1} == 0, not zeros(?,?)
            %elseif numel(x) == 0
            %    x = 0;
            %end
            u = [u x];
            idx = 1:L0;
            listprob = listprob + Lp(:,f{1}+1); % Or better avoid modifing it?
        end
    else
        % Compute soft mapping back one stage
        u1est = cnop(y(:,1:2:end),y(:,2:2:end));
        % R_N^T maps u1est to top polar code
        [uhat1,u1hardprev,idx1,prob1] = polar_decode_list(u1est,f(1:(N/2)),L,listprob,u);
        % Using u1est and x1hard, we can estimate u2
        %u2est = vnop(cnop(u1hardprev,y(1:2:end)),y(2:2:end));
        %u2est = vnop(cnop(u1hardprev,y(idx1,1:2:end)),y(idx1,2:2:end));
        u2est = vnop(y(idx1,1:2:end), y(idx1,2:2:end), u1hardprev);
        % R_N^T maps u2est to bottom polar code
        [uhat2,u2hardprev,idx2,prob2] = polar_decode_list(u2est,f((N/2+1):end),L,prob1,uhat1);
        % Tunnel u decisions back up. Compute and interleave x1,x2 hard decisions
        %u = [uhat1(idx2,:) uhat2];
        u = uhat2;
        %x = reshape([cnop(u1hardprev(idx2,:),u2hardprev); u2hardprev],1,[]);
        x = zeros(length(idx2), size(u2hardprev,2)*2);
        x(:,1:2:end) = bitxor(u1hardprev(idx2,:),u2hardprev);
        x(:,2:2:end) = u2hardprev;
        listprob = prob2;
        idx = idx1(idx2);
    end
    return
% Check-node operation in P1 domain
% function z=cnop(w1,w2)
%     z = w1.*(1-w2) + w2.*(1-w1);
%     return

function l = cnop(l1,l2)
l = (sign(l1).*sign(l2)).* min(abs(l1),abs(l2));
return

% Bit-node operation in P1 domain
% function z=vnop(w1,w2)
%     z = w1.*w2 ./ (w1.*w2 + (1-w1).*(1-w2));
%     return

function l = vnop(l1, l2, bit)
l = (1 - 2*bit).*l1 + l2;
return

function LogP = llr2logp(llr)
denom=log1p(exp(llr));
denom(llr>300)=llr(llr>300);
logp0 = llr-denom;
logp1 = -denom;
LogP=[logp0 logp1];
return

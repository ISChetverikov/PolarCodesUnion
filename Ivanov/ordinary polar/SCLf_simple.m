function [u,x,idx,listprob] = SCLf_simple(y,f,L,listprob,u,inv)
if nargin < 5, u = []; end
if nargin < 4, listprob = 0; end
if nargin < 3, L = 1; end
if ~iscell(f), f = num2cell(f); end
N = size(y,2);
if (N==1)
    Lp=[(y<0).*y,-(y>=0).*y];  %simple approximation
    L0 = size(y,1);
    if (f{1}==1/2)
        % Make list of all possible decisions (double the list size)
        [listprob,j] = sort([listprob + Lp(:,2); listprob + Lp(:,1)],1,'descend');
        % Trim list if it is too large
        if length(j) > L
            if inv~=1
                listprob = listprob(1:L);
                j = j(1:L);
            else
                listprob = listprob(L+1:end);
                j = j(L+1:end);
            end
        end
        % the first half of the list was due to x=1
        x = j <= L0;
        idx = mod(j-1,L0)+1;
        u = [u(idx,:) x];
        % determine which path each decision belongs to
    else
        x = repmat(f{1},size(y));
        u = [u x];
        idx = 1:L0;
        listprob = listprob + Lp(:,f{1}+1); % Or better avoid modifing it?
    end
    
elseif all(cell2mat(f)==1/2)
    [u,x,idx,listprob]=decRate1(y,L,listprob,u,inv); %dec Rate1
elseif all(cell2mat(f)==0)
    [u,x,idx,listprob]=decRate0(y,listprob,u); %dec Rate0
elseif f{end} == 1/2 && all(cell2mat(f(1:end-1)) == 0)
    if isempty(u)
        u=zeros(1,0);
    end
    [u,x,idx,listprob]=decRepCode(y,L,listprob,u,inv); %dec Repetition code
elseif f{1} == 0 && all(cell2mat(f(2:end)) == 1/2)
    [u,x,idx,listprob]=decSPCCode(y,L,listprob,u,inv); %dec SPC code
else
    % Compute soft mapping back one stage
    u1est = cnop(y(:,1:2:end),y(:,2:2:end));
    % R_N^T maps u1est to top polar code
    [uhat1,u1hardprev,idx1,prob1] = SCLf_simple(u1est,f(1:(N/2)),L,listprob,u,inv);
    % Using u1est and x1hard, we can estimate u2
    %u2est = vnop(cnop(u1hardprev,y(1:2:end)),y(2:2:end));
    %u2est = vnop(cnop(u1hardprev,y(idx1,1:2:end)),y(idx1,2:2:end));
    u2est = vnop(y(idx1,1:2:end), y(idx1,2:2:end), u1hardprev);
    % R_N^T maps u2est to bottom polar code
    [uhat2,u2hardprev,idx2,prob2] = SCLf_simple(u2est,f((N/2+1):end),L,prob1,uhat1,inv-N/2);
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
end

% Check-node operation in P1 domain
% function z=cnop(w1,w2)
%     z = w1.*(1-w2) + w2.*(1-w1);
%     return

function l = cnop(l1,l2)
l = (sign(l1).*sign(l2)).* min(abs(l1),abs(l2));
end

% Bit-node operation in P1 domain
% function z=vnop(w1,w2)
%     z = w1.*w2 ./ (w1.*w2 + (1-w1).*(1-w2));
%     return

function l = vnop(l1, l2, bit)
l = (1 - 2*bit).*l1 + l2;
end

function LogP = llr2logp(llr)
denom=log1p(exp(llr));
denom(llr>300)=llr(llr>300);
logp0 = llr-denom;
logp1 = -denom;
LogP=[logp0 logp1];
end

function [u,x,idx,listprob]=decRate1(y,L,listprob,u,inv)
y1=y;
u_hat=kron(y<0,ones(4,1));
logp=-log1p(exp(-abs(y)));
listprob=kron(listprob+sum(logp,2),ones(4,1));

for i = 1:size(y,1)
    [min_val,min_id]=sort(abs(y1(i,:)));
    s=(i-1)*4;
    u_hat([2+s,4+s],min_id(1))=~u_hat([2+s,4+s],min_id(1));
    u_hat([3+s,4+s],min_id(2))=~u_hat([3+s,4+s],min_id(2));
    listprob([2+s,4+s])=listprob([2+s,4+s])-min_val(1);
    listprob([3+s,4+s])=listprob([3+s,4+s])-min_val(2);
end
[listprob,id]=sort(listprob,'descend');
index_list=min(L,size(listprob,1));
if (inv>0) && (inv<=size(u_hat,2))
    listprob=listprob(index_list+1:end);
    id=id(index_list+1:end);
else
    listprob=listprob(1:index_list);
    id=id(1:index_list);
end
x=u_hat(id,:);
idx=floor((id-1)/4)+1;
u=[u(idx,:) polar_transform_mat(x)];
end

function [u,x,idx,listprob]=decRate0(y,listprob,u)
N=size(y,2);
LogP = llr2logp(y);
listprob=listprob+sum(LogP(:,1:N),2);
idx=1:length(listprob);
x=zeros(size(y));
if size(u,2)==0
    u = zeros(size(y));
else
    u=[u(idx,:) x];
end
end

function [u,x,idx,listprob]=decRepCode(y,L,listprob,u,inv)
N=size(y,2);
u_hat=[zeros(size(y)); ones(size(y))];
LogP = llr2logp(y);
listprob=[listprob+sum(LogP(:,1:N),2); listprob+sum(LogP(:,N+1:end),2)];
[listprob,id]=sort(listprob,'descend');
index_list=min(L,size(listprob,1));
if (inv>0) && (inv<=size(u_hat,2))
    listprob=listprob(index_list+1:end);
    id=id(index_list+1:end);
else
    listprob=listprob(1:index_list);
    id=id(1:index_list);
end
x=u_hat(id,:);
idx = mod(id-1, size(y,1)) + 1;
u=[u(idx,:) polar_transform_mat(x)];
end

function [u,x,idx,listprob]=decSPCCode(y,L,listprob,u,inv)
y1=y;
f=0;
u_hat=kron(y<0,ones(8,1));
logp=-log1p(exp(-abs(y)));
listprob=kron(listprob+sum(logp,2),ones(8,1));
for i = 1:size(y,1)
    q=mod(sum(y(i,:)<0)+f,2);
    [min_val,min_id]=sort(abs(y1(i,:)));
    s=(i-1)*8;
    if q == 1
        u_hat(s+1:s+8,min_id(1))=~u_hat(s+1:s+8,min_id(1));
        listprob(s+1:s+8)=listprob(s+1:s+8)-min_val(1);
        min_val(1)=-min_val(1);
    end
    V=nchoosek(1:4,2);
    for j=1:size(V,1)
        u_hat(s+j+1,min_id(V(j,:)))=~u_hat(s+j+1,min_id(V(j,:)));
        listprob(s+j+1)=listprob(s+j+1)-sum(min_val(V(j,:)));
    end
    u_hat(s+8,min_id(1:4))=~u_hat(s+8,min_id(1:4));
    listprob(s+8)=listprob(s+8)-sum(min_val(1:4));
end
[listprob,id]=sort(listprob,'descend');
index_list=min(L,size(listprob,1));

if (inv>0) && (inv<=size(u_hat,2))
   listprob=listprob(index_list+1:end);
    id=id(index_list+1:end); 
else
    listprob=listprob(1:index_list);
    id=id(1:index_list);
end

x=u_hat(id,:);
idx=floor((id-1)/8)+1;
u=[u(idx,:) polar_transform_mat(x)];
end



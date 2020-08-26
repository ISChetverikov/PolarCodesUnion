function [u,list_u,x,idx,listprob,prev] = SCLf_new(y,f,L,listprob,u,inv,known_u,known_listprob)
if nargin < 8, known_listprob = []; end
if nargin < 7, known_u = []; end
if nargin < 6, inv = -1; end
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
        idx = (1:L0)';
        listprob = listprob + Lp(:,f{1}+1); % Or better avoid modifing it?
    end
    prev=[idx; repmat(-1,L-length(idx),1)];
    list_u=u(:,end);
    %      elseif all(cell2mat(f)==1/2)
    %             [u,x,idx,listprob,listprob_alpha]=decRate1(y,L,listprob,listprob_alpha,u,inv,alpha); %dec Rate1
    %      elseif all(cell2mat(f)==0)
    %              [u,x,idx,listprob,listprob_alpha]=decRate0(y,listprob,listprob_alpha,u,inv,alpha); %dec Rate0
    %      elseif f{end} == 1/2 && all(cell2mat(f(1:end-1)) == 0)
    %               [u,x,idx,listprob,listprob_alpha]=decRepCode(y,L,listprob,listprob_alpha,u,inv,alpha); %dec Repetition code
    %      elseif f{1} == 0 && all(cell2mat(f(2:end)) == 1/2)
    %               [u,x,idx,listprob,listprob_alpha]=decSPCCode(y,L,listprob,listprob_alpha,u,inv,alpha); %dec SPC code
else
    if size(known_u,2)>=N/2
        uhat1=[u known_u(:,1:N/2)];
        list_u1=known_u(:,1:N/2);
        u1hardprev=polar_transform_mat(known_u(:,1:N/2));
        if size(y,1)>1
            idx1=1:L;
        else
            idx1=ones(L,1);
        end
        prob1=known_listprob(:,1:N/2);
        prev1=zeros(size(known_u(:,1:N/2)));
    else
        u1est = cnop(y(:,1:2:end),y(:,2:2:end));
        [uhat1,list_u1,u1hardprev,idx1,prob1,prev1] = SCLf_new(u1est,f(1:(N/2)),L,listprob,u,inv,known_u,known_listprob); 
    end
    u2est = vnop(y(idx1,1:2:end), y(idx1,2:2:end), u1hardprev);
    [uhat2,list_u2,u2hardprev,idx2,prob2,prev2] = SCLf_new(u2est,f((N/2+1):end),L,prob1(:,end),uhat1,inv-N/2,known_u(:,N/2+1:end),known_listprob(:,N/2+1:end));
    u = uhat2;
    x = zeros(length(idx2), size(u2hardprev,2)*2);
    x(:,1:2:end) = bitxor(u1hardprev(idx2,:),u2hardprev);
    x(:,2:2:end) = u2hardprev;
    listprob = [[prob1; repmat(-1,size(prob2,1)-size(prob1,1),size(prob1,2))] prob2];
    list_u = [[list_u1; repmat(-1,size(list_u2,1)-size(list_u1,1),size(list_u1,2))] list_u2]; 
    idx = idx1(idx2);
    prev=[prev1 prev2];
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

function [u,x,idx,listprob,listprob_alpha,Malpha]=decRate1(y,L,listprob,listprob_alpha,u,inv,alpha)
u_hat=kron(y<0,ones(4,1));
if (inv>0) && (inv<=size(u_hat,2))
    u_hat(inv,:)=~u_hat(inv,:);
end
logp=-log1p(exp(-abs(y)));
logpa=-log1p(exp(-alpha*abs(y)));
listprob=kron(listprob+sum(logp,2),ones(4,1));
listprob_alpha=kron(listprob_alpha+sum(logpa,2),ones(4,1));
for i = 1:size(y,1)
    [min_val,min_id]=sort(abs(y(i,:)));
    s=(i-1)*4;
    u_hat([2+s,4+s],min_id(1))=~u_hat([2+s,4+s],min_id(1));
    u_hat([3+s,4+s],min_id(2))=~u_hat([3+s,4+s],min_id(2));
    listprob([2+s,4+s])=listprob([2+s,4+s])-min_val(1);
    listprob([3+s,4+s])=listprob([3+s,4+s])-min_val(2);
end
[listprob,id]=sort(listprob,'descend');
index_list=min(L,size(listprob,1));
listprob=listprob(1:index_list);
listprob_alpha=listprob_alpha(1:index_list);
Malpha=cummax(listprob_alpha(id(1:index_list)));
Malpha=Malpha(end-size(u_hat,2):end);
id=id(1:index_list);
x=u_hat(id,:);
idx=floor((id-1)/4)+1;
u=[u(idx,:) polar_transform_mat(x)];
end

function [u,x,idx,listprob,listprob_alpha]=decRate0(y,listprob,listprob_alpha,u,inv,alpha)
N=size(y,2);
LogP = llr2logp(y);
LogPa = llr2logp(alpha*y);
listprob=listprob+sum(LogP(:,1:N),2);
idx=1:length(listprob);
x=zeros(size(y));
if size(u,2)==0
    u = zeros(size(y));
else
    u=[u(idx,:) x];
end
end

function [u,x,idx,listprob]=decRepCode(y,L,listprob,u)
N=size(y,2);
u_hat=[zeros(size(y)); ones(size(y))];
LogP = llr2logp(y);
listprob=[listprob+sum(LogP(:,1:N),2); listprob+sum(LogP(:,N+1:end),2)];
[listprob,id]=sort(listprob,'descend');
index_list=min(L,size(listprob,1));
listprob=listprob(1:index_list);
id=id(1:index_list);
x=u_hat(id,:);
idx = mod(id-1, size(y,1)) + 1;
u=[u(idx,:) polar_transform_mat(x)];
end

function [u,x,idx,listprob]=decSPCCode(y,L,listprob,u)
u_hat=kron(y<0,ones(8,1));
logp=-log1p(exp(-abs(y)));
listprob=kron(listprob+sum(logp,2),ones(8,1));
for i = 1:size(y,1)
    q=mod(sum(y(i,:)<0),2);
    [min_val,min_id]=sort(abs(y(i,:)));
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
listprob=listprob(1:index_list);
id=id(1:index_list);
x=u_hat(id,:);
idx=floor((id-1)/8)+1;
u=[u(idx,:) polar_transform_mat(x)];
end



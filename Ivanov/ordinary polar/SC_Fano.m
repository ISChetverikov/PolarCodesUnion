function u = SC_Fano(y,frozen_bits_indicator,delta,p_e)
A=find(frozen_bits_indicator==0); %information set
frozen_bits_indicator_ext=frozen_bits_indicator;  % for positions previous decoded
N=length(y); %length of codeword
K=length(A); %length of inf
i=1;
j=0;
B=0;
T=0;
u=zeros(1,N); %u_hat - constructed codeword
beta=zeros(1,K);
gamma=zeros(1,K);
path_SC=zeros(1,N); %path metrics
frozen_values=zeros(1,N); %initial frozen values are zeros
while i~=(N+1) %main loop
    if (i-1>0) % we have at least one decision
        if ismember(i-1,A) %information bit
            path_SC(i-1)=L(u(i-1)+1)-log(1-p_e(i-1)); %recalculate path
        else
            path_SC(i-1)=L(1)-log(1-p_e(i-1)); % really frozen bit
        end
    end
    if ismember(i,A) %i in A
        frozen_bits_indicator_ext=frozen_bits_indicator;
        frozen_bits_indicator_ext(1:i-1)=1; %extend frozen positions to previously decoded bits
        frozen_values(1:i-1)=u(1:i-1); %extend frozen positions to previously decoded bits
        [~,~,llr] = polar_decode_llr_stable(y, frozen_bits_indicator_ext,frozen_values,[],-1); % SC decoding to obtain log P(u_i |y, u(1..i-1))
        L=llr2logp(llr(i)); %transform LLR to logP
        %disp(['i=',num2str(i)]) %just for representation
        %disp(['path_i=',num2str(sum(path_SC(1:i-1)))]) %just for representation
        %disp(['T=',num2str(T)]) %just for representation
        mi0=sum(path_SC(1:i-1))+L(1)-log(1-p_e(i));
        mi1=sum(path_SC(1:i-1))+L(2)-log(1-p_e(i));
        if max(mi0,mi1)>T
            if B == 0
                [bb,ind]=max([mi0, mi1]);
                u(i)=ind-1;
                beta(j+1)=bb;
                gamma(j+1)=0;
                if j == 0       
                    mu=0;
                else
                    mu=beta(j);
                end
                if mu<T+delta
                    T=upd(T,delta,beta(j+1));
                end
                i=i+1;
                j=j+1;
            else
                if min([mi0, mi1])>T
                    [mm,idm]=min([mi0, mi1]);
                    u(i)=idm-1;
                    beta(j+1)=mm;
                    gamma(j+1)=1;
                    i=i+1;
                    j=j+1;
                    B=0;
                else
                    if j == 0
                        T=T-delta;
                        B=0;
                    else
                        [T,j,B]=backward(beta,j,T,gamma,delta);
                        i=A(j+1);
                    end
                end
            end
        else
            if j == 0
                T=T-delta;
            else
                [T,j,B]=backward(beta,j,T,gamma,delta);
                i=A(j+1);
            end
        end
    else %i is really frozen
        u(i)=0;
        [~,~,llr_fr] = polar_decode_llr_stable(y, frozen_bits_indicator_ext,frozen_values,[],-1);  % SC decoding to obtain log P(0       |y, u(1..i-1))
        L=llr2logp(llr_fr(i)); %get logP(0)
        i=i+1;
    end
end
end


function T=upd(T,dlt,bb)
while T+dlt < bb
    T=T+dlt;
end
end

function [T,j,B]=backward(beta,j,T,gamma,delta)
while 1
    if j == 1
        mu=0;
    elseif (j>=2)
        mu=beta(j-1);
    end
    if mu>=T
        j=j-1;
        if (gamma(j+1) == 0)    || (j == 0) %?????
            B=1;
            return
        end
    else
        T=T-delta;
        B=0;
        return
    end
end
end

function LogP = llr2logp(llr)
if llr>308
    denom=llr;
elseif llr<-308
    denom=0;
else
    denom=log1p(exp(llr));
end
logp0 = llr-denom;
logp1 = -denom;
LogP=[logp0 logp1];
end


function u = polar_decode_flip_w2(y, frozen_bits_indicator,T1,T21,T22,alpha,crc_poly)
crcdetector = comm.CRCDetector(crc_poly);
N=length(y);
[u_1, ~, out_llr] = polar_decode_llr_stable(y, frozen_bits_indicator,zeros(1,N),[],-1); %without flipping
[~,err]=crcdetector(u_1(~frozen_bits_indicator)');
L2=[];
if err==0
    u=u_1;
    return
else
    L1=FlipDetermine(out_llr,frozen_bits_indicator,-1,T1,alpha);
    for j=0:T1-1
        [u_hat_k, ~,out_llr]=polar_decode_llr_stable(y, frozen_bits_indicator,zeros(1,N),[],L1(j+1));
        [~,err]=crcdetector(u_hat_k(~frozen_bits_indicator)');
        if err==0
           u=u_hat_k;
           return
        end
        if j<T21
            L2=[L2 FlipDetermine(out_llr,frozen_bits_indicator,L1(j+1),T22,alpha)];
        end
    end
    for i=1:T21
        for j=1:T22
            [u_hat_k, ~,~]=polar_decode_llr_stable(y, frozen_bits_indicator,zeros(1,N),[],[L1(i) L2(j+(i-1)*T22)]);
            [~,err]=crcdetector(u_hat_k(~frozen_bits_indicator)');
            if err==0
                u=u_hat_k;
                return
            end
        end
    end    
end
if err==1
    u=randi(2,1,N)-1;
end
end

function J=FlipDetermine(out_llr,frozen_bits_indicator,k1,T,alpha)
    N=length(frozen_bits_indicator);
    M=zeros(1,N);
    for i=1:N
        if (frozen_bits_indicator(i)==0)&&(i>k1)
              M(i)=Malpha(out_llr(1:i),frozen_bits_indicator,alpha);
              %M(i)=abs(out_llr(i));
        else
            M(i)=Inf;
        end
    end
    %[~,JJ]=sort(M,'descend');
    [~,JJ]=sort(M);
    J=JJ(1:T);    
end

%???????    
function val = Malpha(out_llr,frozen_bits_indicator,alpha)
    S=0;
    for i=1:length(out_llr)-1
        if (frozen_bits_indicator(i)==0)
            S=S+log1p(exp(-alpha*abs(out_llr(i))));
        end
    end
    val=(1/alpha)*log1p(alpha*abs(out_llr(end)))+(1/alpha)*S;
end


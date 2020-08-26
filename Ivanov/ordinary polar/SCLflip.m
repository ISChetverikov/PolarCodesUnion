function [u_hat,Idx_smallest,idx_last] = SCLflip(y,f,L,T,alpha,crc_poly,RR)
crcdetector = comm.CRCDetector(crc_poly);
[decoded_bits,~,~,~,~,Malpha] = SCLf(y,f,L,0,0,[],-1,alpha);
Idx_smallest=[];
idx_last=[];
K=length(find(f==1/2));
u_hat=decoded_bits(1,:);
cands=[];
for l=1:L
    uu=decoded_bits(l,1:end);
    %[~,err]=crcdetector(uu(f>0)'); for matlab2020
    [~,err]=step(crcdetector, uu(f>0)');
    if err == 0
        cands=[cands; uu];
    end
end
if ~isempty(cands)
    u_hat=cands(1,1:end);
    return
else
    if isempty(RR)
        [~,idx]=sort(Malpha,'descend');
        %Idx_smallest=idx(1:min(T,length(find(f==1/2))));
        [a, ind]=intersect(idx,find(f==1/2));
        uu=[ind'; a];
        idx=sortrows(uu')';
        Idx_smallest=idx(2,1:min([T,length(a),length(find(f==1/2))]));
    else
        [a, ind]=intersect(RR(2,:),find(f==1/2));
        uu=[ind'; a];
        idx=sortrows(uu')';
        Idx_smallest=idx(2,1:min([T,length(a),length(find(f==1/2))]));
    end    
    for i=1:min([T,K,length(Idx_smallest)])
        decoded_bits = SCLf(y,f,L,0,0,[],Idx_smallest(i),alpha);
        for l=1:L
            uu=decoded_bits(l,1:end);
            %[~,err]=crcdetector(uu(f>0)'); ; for matlab2020
            [~,err]=step(crcdetector, uu(f>0)');
            if err == 0
                u_hat=uu;
                idx_last=Idx_smallest(i);
                return;               
            end
        end      
    end
end
end


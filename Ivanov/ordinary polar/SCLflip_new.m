function [u_hat,Idx_smallest,idx_last] = SCLflip_new(y,f,L,T,crc_poly,RR)
crcdetector = comm.CRCDetector(crc_poly);
[decoded_bits,list_u,~,~,listprob,prev] = SCLf_new(y,f,L,0,[],-1);
Idx_smallest=[];
idx_last=[];
K=length(find(f==0.5));
u_hat=decoded_bits(1,:);
cands=[];
for l=1:L
    uu=decoded_bits(l,1:end);
    [~,err]=step(crcdetector, uu(f>0)');
    if err == 0
        cands=[cands; uu];
    end
end
if ~isempty(cands)
    u_hat=cands(1,1:end);
    return
else
    [a, ind]=intersect(RR(2,:),find(f==1/2));
    uu=[ind'; a];
    idx=sortrows(uu')';
    Idx_smallest=idx(2,1:min([T,length(a),length(find(f==1/2))]));
    for i=1:min([T,K,length(Idx_smallest)])
        known_u=recover_from_list(list_u(:,1:Idx_smallest(i)-1),prev(:,1:Idx_smallest(i)-1));
        known_listprob=listprob(:,1:Idx_smallest(i)-1);
        decoded_bits = SCLf_new(y,f,L,0,[],Idx_smallest(i),known_u,known_listprob);
        for l=1:L
            uu=decoded_bits(l,1:end);
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
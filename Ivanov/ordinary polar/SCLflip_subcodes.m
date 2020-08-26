function [u_hat,Idx_smallest,idx_last] = SCLflip_subcodes(y,f,L,T,Tinner,crc_poly,RR,R1,SPC,Rep)
crcdetector = comm.CRCDetector(crc_poly);
[decoded_bits] = SCLf_simple(y,f,L,0,[],-1);

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
    %new Huawei idea
    if Tinner>0
        MasR1 = find_flip(R1,Idx_smallest);
        MasSPC = find_flip(SPC,Idx_smallest);
        MasREP = find_flip(Rep,Idx_smallest);

        MrelR1 = add_reliability(MasR1,y,f);
        MrelSPC = add_reliability(MasSPC,y,f);
        MrelREP = add_reliability(MasREP,y,f);

        ALL=[MrelR1; MrelSPC; MrelREP];
        ALL=ALL(:,3:4)';
        [r,idx]=sort(ALL(2,:));
        AA=[ALL(1,idx); r];
        Idx_smallest=AA(1,1:min([size(AA,2),T,Tinner,length(a),length(find(f==1/2))]));
    end    
    for i=1:min([T,K,length(Idx_smallest)])
        decoded_bits = SCLf_simple(y,f,L,0,[],Idx_smallest(i));
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


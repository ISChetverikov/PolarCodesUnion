function u_hat = SCL_SP(y,f,L,crc_poly,CS)
crcdetector = comm.CRCDetector(crc_poly);
[decoded_bits] = SCLf(y,f,L,0,0,[],-1,0.375);
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
    CS=intersect(CS,find(f==1/2));
    T=length(CS);
    for i=1:min([T,K])
        decoded_bits = SCLf(y,f,L,0,0,[],CS(i),0.375);
        for l=1:L
            uu=decoded_bits(l,1:end);
            [~,err]=step(crcdetector, uu(f>0)');
            if err == 0
                u_hat=uu;
                return;               
            end
        end      
    end
end
end


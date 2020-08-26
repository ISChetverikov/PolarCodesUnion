function [cands,L] = adaptive_SCL_decode(y,f,Lmax,crc_poly)
L=1;
frozen_bits_indicator=zeros(1,length(y));
frozen_bits_indicator(f==0)=1;
cands=[];
err=1;
while (L<=Lmax) && (err==1)
    crcdetector = comm.CRCDetector(crc_poly);
    decoded_bits = polar_decode_list(y,f,L);
    for l=1:L
        uu=decoded_bits(l,1:end);
        [~,err]=crcdetector(uu(~frozen_bits_indicator)');
        if err == 0
            cands=[cands; uu];
            return;
        end
    end
    if isempty(cands)
        L=L*2;
    end
end
if err == 1
    cands = randi(2,1,length(y))-1;
end   
end


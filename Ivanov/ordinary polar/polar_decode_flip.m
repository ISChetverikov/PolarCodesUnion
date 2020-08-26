function u = polar_decode_flip(y, frozen_bits_indicator,T,crc_poly)
crcdetector = comm.CRCDetector(crc_poly);
[u_1, ~, out_llr] = polar_decode(y, frozen_bits_indicator,[],length(y)+1); %without flipping
[~,err]=crcdetector(u_1(~frozen_bits_indicator)');
if err==0
    u=u_1;
else
    I=find(frozen_bits_indicator==0); %information set positions
    if (T>1) && (err==1) %if err==1?
        B=[out_llr(I); I]; %llrs in these positions
        [~, inx]=sort(abs(B(1,:))); %sort by LLR
        SortedMat=B(:,inx); 
        U=SortedMat(2,1:min(T,length(SortedMat))); %form set of T least reliable channels
        for j=1:T
            k=U(j);
            [u_hat_k, ~,~]=polar_decode(y, frozen_bits_indicator,[],k);
            [~,err]=crcdetector(u_hat_k(~frozen_bits_indicator)');
            if err==0
                %disp(sprintf('succ on %d iter',j));
                u=u_hat_k;
                return
            end            
        end
    %disp('error decoding');
    u=y<0;
    return
    end
end
end


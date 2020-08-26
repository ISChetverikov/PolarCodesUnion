function u_hat = polar_decode_chase(llr, f,chase,crc_poly,Mr)
crcdetector = comm.CRCDetector(crc_poly);
frozen_bits_indicator=~f;
u1 = polar_decode(llr, ~f,[],-1);
[~,err]=step(crcdetector, u1(f>0)');
if err==0
    u_hat=u1;
else
    if chase>0
        llr_cur=llr;
        llr_cur(find(f==0))=Inf;
        [~,idx]=sort(abs(llr_cur));
        a = repmat(llr, 2^chase, 1);
        [cols{1:chase}] = ndgrid([-100 100]);
        a(:, idx(1:chase)) = reshape([cols{:}], [], chase);
        new_LLR=a;
        U=[];
        for i=1:size(new_LLR,1)
            uu = BP_Decoder_LLR((~frozen_bits_indicator)', frozen_bits_indicator', new_LLR(i,:), 200, Mr);
            uu=double(uu);
            [~,err]=step(crcdetector, uu);
            if err == 0
                U=[U; uu];
            end
        end

        if ~isempty(U)
            u_hat=U(1,1:end);
        else
            u_hat=zeros(1,length(llr));
        end
    else
        u_hat=zeros(1,length(llr));
    end
end






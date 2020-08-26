function u = ISD_SCdec(y,frozen_bits_indicator,alpha,crc_poly,Gen)
crcdetector = comm.CRCDetector(crc_poly);
N=length(y);
%[u_1, x_1, out_llr] = polar_decode_llr_stable(y, frozen_bits_indicator,zeros(1,N),[],-1); %without flipping
[info_esti, x_1, llr_out] = BP_Decoder_LLR((~frozen_bits_indicator)', frozen_bits_indicator', y, 50, index_Matrix(N));
%[~,err]=crcdetector(u_1(~frozen_bits_indicator)');
[~,err]=crcdetector(info_esti);
k=size(Gen,1);
if err==0
    %u=u_1(~frozen_bits_indicator);
    u=info_esti';
    return
else
    x_1=x_1';
    llr_out=llr_out';
%        M=zeros(1,N);
%         for j=1:N
%             M(j)=Malpha(out_llr(1:j),frozen_bits_indicator,alpha);
%         end
    %[~,JJ]=sort(-abs(out_llr));
    [~,JJ]=sort(-abs(llr_out));
    %[~,JJ]=sort(M);
    JJk=sort(JJ(1:k));
    if det(gf(Gen(:,JJk)))==1
        u=x_1(JJk)*inv(gf(Gen(:,JJk)));
        return
    else
        [~,jb] = gfrref(gf(Gen(:,JJk)),2);
        rank=length(jb);
        r_def=size(gf(Gen(:,JJk)),2)-rank;
        V=nchoosek(JJ(k+1:N),r_def);
        for l=1:size(V,1)
            Idx=sort([JJk(jb) V(l,:)]);
            Gl=gf(Gen(:,Idx));
            if det(Gl)==1
                u=x_1(Idx)*inv(Gl);
            end
        end
    end
end
end

function val = Malpha(out_llr,frozen_bits_indicator,alpha)
    S=0;
    for i=1:length(out_llr)-1
        if (frozen_bits_indicator(i)==0)
            S=S+log1p(exp(-alpha*abs(out_llr(i))));
        end
    end
    val=(1/alpha)*log1p(alpha*abs(out_llr(end)))+(1/alpha)*S;
end

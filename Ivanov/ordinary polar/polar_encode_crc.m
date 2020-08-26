function [x,u1] = polar_encode_crc(u_inf,crc_poly,frozen_bits_indicator)
crcgenerator = comm.CRCGenerator(crc_poly);
% 'z^16 + z^14 + z + 1' - CRC-16
%u_inner=crcgenerator(u_inf')'; %len(u_inf)=A, length(u_inner)=K=L+A,
%len(CRC)=L 2020? works
u_inner=step(crcgenerator, u_inf')';
N=length(frozen_bits_indicator); %sum(frozen_bits_indicator)=N-K
if sum(frozen_bits_indicator)~=(N-length(u_inner))
    fprintf('error')
    x=-1*ones(1,N);
else
    c=zeros(1,N);
    c(~frozen_bits_indicator)=u_inner;
    x = polar_encode_recursive(c);
    u1=c;
end
end


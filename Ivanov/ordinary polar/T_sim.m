clc;
clear all;
crc_poly='z^6 + z + 1';
crc_poly=[1 0 0 0 0 1 1];
poly_deg=6;
crcdetector = comm.CRCDetector(crc_poly);
crcgenerator = comm.CRCGenerator(crc_poly);
n_array=[7 8 9 10];
R_array=[1/9 1/5 1/3 1/2 2/3];
T_array=[4 8 16 32 64];
L_array=[8 16 32];
Fr = table2array(readtable('rank.csv'));
channels=Fr(:,2)';
M=2;
fe=100;
EbNoArray=-7:0.2:9;
List_size=8;
fer_out=1e-3;
rng('default');
parfor rr = 1:length(R_array)
    R=R_array(rr);
    for nn=1:length(n_array)
        n=n_array(nn);
        N = 2^n;
        K=floor(N*R);
        channels=Fr(:,2)';
        channels=channels(channels<N);
        frozen_indices = channels(1:N-K-poly_deg)+1;
        frozen_bits_indicator=zeros(1,N);
        frozen_bits_indicator(frozen_indices)=1;
        f=frozen_bits_indicator;
        f(frozen_bits_indicator==0)=1/2;
        f(f==1)=0;
        for tt=1:length(T_array)
            T=T_array(tt);
            FER_SCLflip=zeros(1,length(EbNoArray));
            fer_real=1;
            i=0;
            while fer_real>fer_out
                i=i+1;
                EbNo=EbNoArray(i);
                noiseVar = 10.^(-EbNo/10);
                f_sclflip=0;
                Itr=0;
                while f_sclflip<fe
                    Itr=Itr+1;
                    u = randi([0, 1], 1, K);
                    [x,u1] = polar_encode_crc(u,crc_poly,frozen_bits_indicator); %encode with CRC
                    txSig = qammod(x,M,'InputType','bit','UnitAveragePower',true);
                    rxSig = awgn(txSig,EbNo,'measured');
                    rxDataSoft = qamdemod(rxSig,M,'OutputType','approxllr', ...
                        'UnitAveragePower',true,'NoiseVariance',noiseVar);
                    decoded_bits8 = SCLflip(rxDataSoft,f,List_size,T,0.375,crc_poly);
                    
                    if ~isequal(decoded_bits8, u1)
                        f_sclflip=f_sclflip+1;
                    end
                    FER_SCLflip(i)=f_sclflip/Itr;
                end
                fer_real=FER_SCLflip(i);
                disp([ 'Block Error Rate SCL flip: ' num2str(FER_SCLflip(i)) ...
                    ', at SNR = ' num2str(EbNo) ' dB. T=' num2str(T) ...
                    ', N =' num2str(N) ...
                    ',R = ' num2str(R)]);
            end
            save(sprintf('SCL_flip_T=%d_N=%d_R=%f.mat', T,N,R));
        end
        for l=1:length(L_array)
            FER_SCL=zeros(1,length(EbNoArray));
            L=L_array(l);
            fer_real=1;
            i=0;
            while fer_real>fer_out
                i=i+1;
                EbNo=EbNoArray(i);
                noiseVar = 10.^(-EbNo/10);
                f_scl=0;
                Itr1=0;
                while f_scl<fe
                    Itr1=Itr1+1;
                    u = randi([0, 1], 1, K);
                    [x,u1] = polar_encode_crc(u,crc_poly,frozen_bits_indicator); %encode with CRC
                    txSig = qammod(x,M,'InputType','bit','UnitAveragePower',true);
                    rxSig = awgn(txSig,EbNo,'measured');
                    rxDataSoft = qamdemod(rxSig,M,'OutputType','approxllr', ...
                        'UnitAveragePower',true,'NoiseVariance',noiseVar);
                    decoded_bits1 = polar_decode_list(rxDataSoft,f,L);
                    
                    if ~isequal(decoded_bits1(1,:), u1)
                        f_scl=f_scl+1;
                    end
                    FER_SCL(i)=f_scl/Itr1;
                end
                fer_real=FER_SCL(i);
                disp([ 'Block Error Rate SCL: ' num2str(FER_SCL(i)) ...
                    ', at SNR = ' num2str(EbNo) ' dB. L=' num2str(L) ...
                    ', N =' num2str(N) ...
                    ',R = ' num2str(R)]);
            end
            save(sprintf('SCL_L=%d_N=%d_R=%f.mat', L,N,R));
        end    
    end
end




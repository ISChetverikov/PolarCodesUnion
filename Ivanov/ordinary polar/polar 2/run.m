%% Configuration
p = 5;
K = 64;
R = 1/4;
do_crc = 1;
list_size = 32;
crc1_idx = 10; %0x1A-10, 0x23-25, 0x9-11
%crc1_polynom = hex2dec('9'); %4
crc1_polynom = hex2dec('1A'); %5
%crc1_polynom = hex2dec('23'); %6
%crc2_polynom = hex2dec('1d8ae'); %17
crc2_polynom = hex2dec('212d'); %14
%crc2_polynom = hex2dec('573a'); %15
%crc2_polynom = hex2dec('9eb2'); %16
% crc1_polynom = hex2dec('38');
%crc2_polynom = hex2dec('6b99c'); %19
terminate_lately = true;
EsNos = -4i:0.25:0;
random_data = false;
% other = 1; % add random coset



%% Initialization
M = round(K/R);
N = 2^ceil(log2(M));
Fr = table2array(readtable('rank.csv'));
channels=Fr(:,2)';
channels=channels(channels<N);

if crc1_idx > 0
    crc1_len = numel(de2bi(crc1_polynom));
else
    crc1_len = 0;
end
if do_crc
    crc2_len = numel(de2bi(crc2_polynom));
    frozen_indices = channels(1:N-K-crc1_len-crc2_len)+1;
    frozen_bits_indicator=zeros(1,N);
    frozen_bits_indicator(frozen_indices)=1;
    ff=frozen_bits_indicator;
    ff(frozen_bits_indicator==0)=1/2;
    ff(ff==1)=0;
    f=char(zeros(1,N));
    f(ff==0.5)='i';
    f(ff==0)='f';
else
    crc2_len = 0;
    frozen_indices = channels(1:N-K-crc1_len)+1;
    frozen_bits_indicator=zeros(1,N);
    frozen_bits_indicator(frozen_indices)=1;
    ff=frozen_bits_indicator;
    ff(frozen_bits_indicator==0)=1/2;
    ff(ff==1)=0;
    f=char(zeros(1,N));
    f(ff==0.5)='i';
    f(ff==0)='f';
end
%f(f=='p') = 'f';
if crc1_idx > 0
    a = find(f=='i');
    crc1_pos = a(crc1_idx + 1 : crc1_idx+crc1_len);
    K1 = crc1_idx;% - crc1_len;
    %f(crc1_pos) = 'c';
else
    K1 = 1;
    crc1_pos = 0;
end
fprintf('[%d,%d] N=%d R=%.2f crc-%d + crc-%d\n', M, K, N, K/M, crc1_len, numel(de2bi(crc2_polynom)));
Errs = zeros(numel(EsNos),5);
Early = zeros(numel(EsNos),1);
TermPos = zeros(numel(EsNos),N);
fname = [char(datetime) '.mat'];
crc_mask = find(f=='i' | f=='c');
mode = {'signal', 'noise'};
%crc_mask(crc1_idx - crc_len+1:crc1_idx) = 0;
%crc_mask = crc_mask(K1+1:end);

%% Generator matrix
G = false(K,N);
G2 = false(K,K+crc1_len+crc2_len);
gen1 = crc.generator('Polynomial', [1 de2bi(crc1_polynom)], 'InitialState', 0); % CRC-6
gen2 = crc.generator('Polynomial', [1 de2bi(crc2_polynom)], 'InitialState', 0); % CRC-19
for i=1:K
    I = zeros(1,K);
    I(i) = 1;
    if crc1_idx > 0
        I = [gen1.generate(I(1:K1)')' I(K1+1:end)];
    end
    sc = zeros(1,N);
    if do_crc
        %sc(f == 'i' | f == 'c') = [I crc19_gen(I(K1+1:end))];
        %sc(f == 'i' | f == 'c') = [I(1:K1) gen2.generate(I(K1+1:end)')'];
        sc(f == 'i' | f == 'c') = gen2.generate(I')';
        G2(i,:) = gen2.generate(I')';
        assert(crc_ok(crc2_polynom, sc(crc_mask)));
    else
        sc(f == 'i' | f == 'c') = I;
    end
    for j=find(f=='p')
        t = mod(j-1,p)+1;
        sc(j) = mod(sum(sc(t:p:j-1)),2);
    end
    G(i,:) = polar_transform(sc);
end
G = G(:,1:M);
final_filename = sprintf('%s [%d,%d] crc%d+crc%d', mode{random_data+1}, M,K, crc1_len, numel(de2bi(crc2_polynom)));
display(final_filename);

% G2_s = find(f=='p' | f=='f');
% G2 = false(length(G2_s),M);
% for i=1:length(G2_s)
%     sc = false(1,M);
%     sc(G2_s(i)) = true;
%     G2(i,:) = polar_transform(sc);
% end
crc1_end = crc1_pos(end);

%% Simulation
for i=1:numel(EsNos)
    EsNo = EsNos(i);
    noise_var = 1/(exp(log(10) * EsNo/10));
    time = tic;
    dec_time = 0;
    e = Errs(i,:);
    early = 0;
    term_pos = TermPos(i,:);
    while e(3) + e(4) < 100 || (random_data && e(5) < 1e6)
        parfor sdf=1:10000
            I = randi([0 1], 1, K);
            if ~random_data
                cw = 1 - 2*mod(I*G,2);
            else
                cw = 1 - 2*randi([0 1], 1, M);
            end
            noise = randn(size(cw))/sqrt(2);
            rx = cw + sqrt(noise_var) * noise;
            llr = 4/noise_var * real(rx);
            llr = [llr zeros(1,N-M)];
            y = 1 ./ (1+exp(llr));
            t = tic;
            [dec, prob] = pcscl(-llr,f,p,list_size, crc1_end, crc1_polynom,[[]],terminate_lately);
            dec_time = dec_time + toc(t);
            if size(prob,2) < N+1
                e = e + [0, 0, 0, 1, 1]; % reject
                early = early + 1;
                term_pos = term_pos + accumarray(size(dec,2),1,[N 1])';
                continue;
            end
            assert(size(dec,2) == N && size(prob, 2) == N+1);
            if do_crc
                crcerr = crc_ok(crc2_polynom,dec(:,crc_mask));
                dec_i = dec(crcerr,f=='i' | f=='c');
            else
                [~,b] = max(prob(:,end));
                dec_i = dec(b,f=='i' | f=='c');
                crcerr = true(1,size(dec,1));
            end
            if size(dec_i,1) < 1
                e = e + [0, 0, 0, 1, 1]; % reject
                fail = 1;
            else
                [~,b] = max(prob(crcerr,end));
                dec_i = [dec_i(b,1:K1) dec_i(b,K1+crc1_len+1:K+crc1_len)];
                fail = ~isequal(dec_i,I);
                e = e + [0, 0, fail, 0, 1];
            end
        end
        Early(i) = early;
        Errs(i,:) = e;
        assert(all(~term_pos(f~='i')));
        TermPos(i,:) = term_pos;
        save(fname);
        fprintf('\tSNR %.2fdB,\tFER: early %3d + error %3d + full %3d/%6d = %.0e\tFAR: %.0e\ttime: %.1f\n', ...
            EsNo, early, e(3), e(4)-early, e(5), (e(3)+e(4))/e(5), e(3)/e(5), toc(time));
    end
    fprintf('SNR %.2fdB,\tFER: early %3d + error %3d + full %3d/%6d = %.0e\tFAR: %.0e\ttime: %.1f, %.1f Kbps\n', ...
        EsNo, early, e(3), e(4)-early, e(5), (e(3)+e(4))/e(5), e(3)/e(5), toc(time), e(5)*K/dec_time/1e3);
    saved_code = sum(term_pos .* (N-1:-1:0)) / e(5);
    saved_info = sum(term_pos(f=='i') .* (nnz(f=='i')-1:-1:0)) / e(5);
    fprintf('N,M,K,rate,crc1-end,snr,early,term-pos,saved-code,saved-info,err,fail,frames,crc1-len,crc2-len,zeroseq\n');
    fprintf('%d,%d,%d,%s,%d,%.1f,%d,%.2f,%.3f,%.3f,%d,%d,%d,%d,%d,%d\n',N,M,K,strtrim(rats(K/M,4)),crc1_end,EsNo,early, sum(term_pos .* (1:N))/Errs(i,4), saved_code, saved_info,Errs(i,3),Errs(i,4),Errs(i,5), crc1_len, crc2_len, find(f=='i',1)-1);
end

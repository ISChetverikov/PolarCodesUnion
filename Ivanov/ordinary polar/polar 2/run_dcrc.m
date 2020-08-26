%% Configuration
p = 5;
K = 64;
R = 1/12;
crc1_len = 3;
crc2_polynom = hex2dec('6b99c');
EsNos = -8;
random_data = true;

%% Initialization
M = round(K/R);
N = 2^ceil(log2(M));
crc2_len = numel(de2bi(crc2_polynom));
f = huawei_f(N,M,K + crc2_len);
[checks_i, perm] = crc_et_checks(crc2_polynom, K, crc1_len);
info_pos = find(f=='i');
crc1_end = info_pos(size(checks_i,2));
checks = zeros(crc1_len, crc1_end);
checks(:,info_pos(1:size(checks_i,2))) = checks_i;
fprintf('[%d,%d] N=%d R=%.2f crc-%d, early termination info %d, code %d\n', M, K, N, K/M, crc2_len, size(checks_i,2), crc1_end);
Errs = zeros(numel(EsNos),5);
Early = zeros(numel(EsNos),1);
fname = [char(datetime) '.mat'];
mode = {'signal', 'noise'};
final_filename = sprintf('%s [%d,%d] crc%d et%d', mode{random_data+1}, M,K,crc2_len, crc1_len);

%% Generator matrix
G = false(K,N);
gen2 = crc.generator([1 de2bi(crc2_polynom)]); % CRC-19
%inv_perm(perm) = 1:length(perm);
clear crc_mask
crc_mask(perm) = find(f=='i' | f=='c');
for i=1:K
    I = zeros(1,K);
    I(i) = 1;
    sc = zeros(1,N);
    I = gen2.generate(I')';
    sc(f == 'i' | f == 'c') = I(perm);
    assert(all(mod(I(perm(1:size(checks_i,2))) * checks_i',2) == 0));
    assert(all(mod(sc(1:crc1_end) * checks',2) == 0));
    assert(crc_ok(crc2_polynom, sc(crc_mask)));
    for j=find(f=='p')
        t = mod(j-1,p)+1;
        sc(j) = mod(sum(sc(t:p:j-1)),2);
    end
    G(i,:) = polar_transform(sc);
end
G = G(:,1:M);

%% Simulation
for i=1:numel(EsNos)
    EsNo = EsNos(i);
    noise_var = 1/(exp(log(10) * EsNo/10));
    time = tic;
    dec_time = 0;
    e = Errs(i,:);
    early = 0;
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
            llr = [4/noise_var * real(rx), zeros(1,N-M)];
            t = tic;
            [dec, prob] = pcscl(-llr,f,p,8, crc1_end, 0, checks, false);
            dec_time = dec_time + toc(t);
            if size(prob,2) == crc1_end
                e = e + [0, 0, 0, 1, 1]; % reject
                early = early + 1;
                continue;
            end
            assert(size(dec,2) == N && size(prob, 2) == N+1);
            crcerr = crc_ok(crc2_polynom,dec(:,crc_mask));
            dec_i = dec(crcerr,crc_mask);
            if size(dec_i,1) < 1
                e = e + [0, 0, 0, 1, 1]; % reject
                fail = 1;
            else
                [~,b] = max(prob(crcerr,end));
                dec_i = dec_i(b,1:K);
                fail = ~isequal(dec_i,I);
                e = e + [0, 0, fail, 0, 1];
            end
        end
        Early(i) = early;
        Errs(i,:) = e;
        save(fname);
        fprintf('\tSNR %.2fdB,\tFER: early %3d + error %3d + full %3d/%6d = %.0e\tFAR: %.0e\ttime: %.1f\n', ...
            EsNo, early, e(3), e(4)-early, e(5), (e(3)+e(4))/e(5), e(3)/e(5), toc(time));
    end
    fprintf('SNR %.2fdB,\tFER: early %3d + error %3d + full %3d/%6d = %.0e\tFAR: %.0e\ttime: %.1f, %.1f Kbps\n', ...
        EsNo, early, e(3), e(4)-early, e(5), (e(3)+e(4))/e(5), e(3)/e(5), toc(time), e(5)*K/dec_time/1e3);
end

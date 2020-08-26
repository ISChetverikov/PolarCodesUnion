function [ f ] = huawei_f( N, M, K )
%HUAWEI_F Summary of this function goes here
%   Detailed explanation goes here
W = de2bi(0:N-1) * 2.^((0:log2(N)-1)'/4);
%assert(N == M);
%P = N-M;
if M<N
    P = bi2de(de2bi(M:N-1,log2(N),'left-msb'),'right-msb')+1;
else
    P = [];
end
alpha = 1.5;
Fp = floor(min(log2(N)*(alpha - (alpha*abs(K/M-0.5))^2), (M-K)/2));
assert(Fp <= M-K);
[~, Q] = sort(W);
Q = [Q(ismember(Q,P)); Q(~ismember(Q,P))];
H = eye(N);
for i=1:N, H(i,:) = polar_transform(H(i,:)); end
RW = sum(H,2);

w_min = min(RW(Q(end-K-Fp+1:end)));
n = nnz(RW(Q(end-K-Fp+1:end)) == w_min);
f1 = ceil((Fp + min(Fp,n))/2);
f2 = floor((Fp - min(Fp,n))/2);
if f1  > n
    f2 = f2 + ceil((f1-n)/2);
    f1 = n;
end
%fprintf('(%d,%d,%d)\n', w_min, f1, f2);
p1 = find(RW(Q) == w_min);
p1 = p1(end-f1+1:end);
p2 = find(RW(Q) == 2*w_min);
p2 = p2(end-f2+1:end);
I = true(N,1);
I(p1) = 0; I(p2) = 0; % disable bit reserved as PC-frozen
%I(P) = 0; % disbale puncutred bits
I = Q(I);
I = I(end-K+1:end);
f = repmat('f',1,N);
f(I) = 'i';
f(and(or(RW == w_min, RW == 2*w_min)', f == 'f')) = 'p';
end


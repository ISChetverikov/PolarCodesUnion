
n=6;
N=64;
G_2 = [1, 0; 1, 1];

G = G_2;

for index = 1:(n-1)
  G = kron(G_2, G);
end
GN=bitrevorder(G);

G=GN(find(f==1/2),:);
H=GN(:,find(f==0));

HH=H';
[srt,id]=sort(sum(HH));
X=[srt; id];
YY=X(:,ismember(X(2,:),find(f==1/2)));
WT=YY(:,ismember(YY(2,:),RR_5G(2,:)));


Fr=Fr';
CH=Fr(:,ismember(Fr(2,:),0:63));

K=corrcoef(CH(1,:),sum(HH))
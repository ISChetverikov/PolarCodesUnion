function CS = generateCS(A,n)
B=-ones(n+1,2^n);
CS=[];
for k=1:2^n
    B(n+1,k)=0;
    if ~ismember(k,A)
        B(n+1,k)=1;
    end
end
for i=n:-1:1
    for j=1:2^(i-1)
        B(i,j)=B(i+1,2*j-1)+B(i+1,2*j);
    end
end
for i=1:(n+1)
    for j=1:2^(i-1)
        x1=j;
        x2=j;
        if B(i,j)==0
            for k=1:(n+1)-i
                x1=2*x1-1;
                x2=2*x2;
                for p=x1:x2
                    B(i+k,p)=-1;
                end
                CS=[CS x1];
            end
        end
    end
end
end


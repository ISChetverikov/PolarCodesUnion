function LogP = llr2logp(llr)
if llr>308
    denom=llr;
elseif llr<-308
    denom=0;
else
    denom=log1p(exp(llr));
end
logp0 = llr-denom;
logp1 = -denom;
LogP=[logp0 logp1];
end

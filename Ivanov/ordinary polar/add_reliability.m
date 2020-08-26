function Mrel = add_reliability(MasT,y,f)
% Mrel=[];
% yy=y;
% yy(f==0)=1e3;
% yy=abs(yy);
% if ~isempty(MasT)
%     for i=1:size(MasT,1)
%         Mrel=[Mrel; [MasT(i,:) min(yy(MasT(i,1):MasT(i,2)))]];
%     end
% end
% end

Mrel=[];
yy=y;
if ~isempty(MasT)
    for i=1:size(MasT,1)
        metric=sum(abs(yy(MasT(i,1):MasT(i,2))).*f(MasT(i,1):MasT(i,2)))/sum(f(MasT(i,1):MasT(i,2)));
        Mrel=[Mrel; [MasT(i,:) metric]];
    end
end
end
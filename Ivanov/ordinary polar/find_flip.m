function MasT = find_flip(Mas,Tvec)
MasT=[];
if ~isempty(Mas)
    for i=1:size(Mas,1)
        if any(ismember(Tvec,Mas(i,1):Mas(i,2)))
            U=Tvec(ismember(Tvec,Mas(i,1):Mas(i,2)));
            MasT=[MasT; [Mas(i,:) U(1)]];
        end
    end
end
end


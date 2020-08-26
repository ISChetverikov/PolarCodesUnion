function [checks,perm] = crc_et_checks(polynom, N, et)
    if ischar(polynom)
        polynom = hex2dec(polynom);
    end
    if ~islogical(polynom)
        polynom = [1 de2bi(polynom)];
    end
    %fprintf('crc-%d\n',length(polynom)-1);
    gen = crc.generator(polynom);
    G = gen.generate(eye(N));
    G = G(N+1:end,:);% * eye(N);
    rows = nchoosek(1:size(G,1),et);
    min_weight = N;
    min_rows = 1:et;
    for i=1:size(rows,1)
        w = nnz(sum(G(rows(i,:),:),1));
        if w < min_weight
            min_weight = w;
            min_rows = rows(i,:);
            checks = G(rows(i,:),:);
            %fprintf('min weight %d, rows: %s\n',w, mat2str(rows(i,:)));
        end
    end
    if nargout > 1
        a = sum(checks,1) > 0;
        other_checks = 1:size(G,1);
        other_checks = other_checks(~ismember(other_checks,min_rows));
        perm = [find(a) (min_rows+N) find(~a) (other_checks+N)];
        checks = [checks(:, a) eye(et)];
    end
end

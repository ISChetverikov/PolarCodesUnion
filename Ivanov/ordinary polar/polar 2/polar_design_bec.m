function f = polar_design_bec(n,e,d)
    % Design polar code of length N=2^n for BEC(e) and target block error rate d
    % Generate virtual channel erasure probabilities
    E = e;
    for i=1:n
        % Interleave updates to keep in polar decoding order
        E = reshape([1-(1-E).*(1-E); E.*E],1,[]);
    end
    % Sort into increasing order and compute cumulative sum
    [SE,order] = sort(E);
    CSE = cumsum(SE);
    % Find good indices
    I = sum(double(CSE<d));
    f = zeros(1,length(E));
    f(order(1:I)) = 1/2;


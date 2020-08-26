function [u, v, out_llr] = polar_decode_terminated(y, u_prev, frozen_bits_indicator,frozen_values,in_llr,stop_pos)
% inv_bit_numb - not in frozen set
N = length(y);
out_llr=in_llr;
if N == 2
    if length(u_prev)<N/2
        L_u_1 = f(y(1), y(2));
        out_llr=[out_llr L_u_1];
        u = zeros(1, 2);
        if frozen_bits_indicator(1) == 1
            u(1) = frozen_values(1);
        elseif L_u_1 >= 0
            u(1) = 0;
        else
            u(1) = 1;
        end
        i=length(out_llr);
        if i == stop_pos
            v=-1;
            return;
        end
    else %not decode left sub-tree
        u(1)=u_prev(1);
        %decode right sub-tree
        L_u_2 = g(y(1), y(2), u(1));
        out_llr=[out_llr L_u_2];
        if frozen_bits_indicator(2) == 1
            u(2) = frozen_values(2);
        elseif L_u_2 >= 0
            u(2) = 0;
        else
            u(2) = 1;
        end
        out_llr(2)=L_u_2;
        i=length(out_llr);
        if i == stop_pos
            v=-1;
            return;
        end
        v = zeros(1, 2);
        v(1) = bitxor(u(1), u(2));
        v(2) = u(2);
    end    
else    
    L_w_odd = zeros(1, N/2);
    for index = 1:(N/2)
        L_w_odd(index) = f(y(2*index-1), y(2*index));
    end
    
    frozen_bits_indicator_1 = frozen_bits_indicator(1:(N/2));
    frozen_values_1 = frozen_values(1:(N/2));

    [u_1, v_1,out_llr] = polar_decode_terminated(L_w_odd, u_prev(1:N/2), frozen_bits_indicator_1,frozen_values_1,out_llr,stop_pos);
    
    L_w_even = zeros(1, N/2);
    for index = 1:(N/2)
        L_w_even(index) = g(y(2*index-1), y(2*index), v_1(index));
    end
    
    frozen_bits_indicator_2 = frozen_bits_indicator((N/2+1):N);
    frozen_values_2 = frozen_values((N/2+1):N);

    [u_2, v_2,out_llr] = polar_decode_terminated(L_w_even, u_prev(N/2+1:end), frozen_bits_indicator_2,frozen_values_2,out_llr,stop_pos);
    
    u = [u_1, u_2];
    
    v = zeros(1, N);
    for index = 1:(N/2)
        v(2*index-1) = bitxor(v_1(index), v_2(index));
        v(2*index) = v_2(index);
    end
    
end

end
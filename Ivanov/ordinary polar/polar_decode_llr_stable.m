function [u,x, out_llr] = polar_decode_llr_stable(llr, F,frozen_values, in_llr,inv_mas)

N = length(llr);
out_llr=in_llr;
if (N==1)
    out_llr=[out_llr llr];
    i=length(out_llr);
    if (F==0)
        u = (llr < 0);
        x = u;
    else
        u = frozen_values(1);
        x = u;
    end
    if ismember(i,inv_mas) %flipping
        u = 1-u;
        x = u;
    end
else
    % Upper part
    llr1 = cnop1(llr(1:2:end), llr(2:2:end));
    [u1, x1,out_llr] = polar_decode_llr_stable(llr1, F(1:(N/2)), frozen_values(1:(N/2)),out_llr,inv_mas);
    
    % Lower part
    llr2 = vnop1(llr(1:2:end), llr(2:2:end), x1);
    [u2, x2,out_llr] = polar_decode_llr_stable(llr2, F((N/2+1):end), frozen_values((N/2+1):end),out_llr,inv_mas);
    
    % Combine
    u = [u1 u2];
    x = reshape([bitxor(x1, x2); x2],1,[]);
end
end

function l = cnop2(l1,l2)
l = (2 * (atanh(tanh(l1/2).*tanh(l2/2))));
end

function l = cnop1(l1,l2)
l = (sign(l1).*sign(l2)).* min(abs(l1),abs(l2));
end


function l = vnop1(l1, l2, bit)
l = ((-1).^bit).*l1 + l2;
end





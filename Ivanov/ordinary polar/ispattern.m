function u = ispattern(f)
    u=[0 0 0 0];
    if all(f==1/2)
        u=[0 1 0 0];
    end
    if all(f==0)
        u=[1 0 0 0];
    end
    if f(end) == 1/2 && all(f(1:end-1) == 0)
        u=[0 0 0 1];
    end
    if f(1) == 0 && all(f(2:end) == 1/2)
        u=[0 0 1 0];
    end     
end

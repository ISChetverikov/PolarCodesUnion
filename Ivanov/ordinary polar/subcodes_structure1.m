function [R0,R1,SPC,Rep] = subcodes_structure1(f)
    N=length(f);
    f_cur=f;
    l_past=N;
    left=0;
    R0=[];
    R1=[];
    SPC=[];
    Rep=[];
    while N>0
        N=length(f_cur);
        pos=find(fix(log2(1:N))-log2(1:N)==0);
        for j=length(pos):-1:1
            if sum(ispattern(f_cur(1:pos(j))))>0
                f_st=f_cur(1:pos(j));
                u=ispattern(f_st);
                if u(1) == 1
                    R0=[R0; left+1 left+length(f_st)];
                    left=left+length(f_st);
                end
                if u(2) == 1
                    R1=[R1; left+1 left+length(f_st)];
                    left=left+length(f_st);
                end
                if u(3) == 1
                    SPC=[SPC; left+1 left+length(f_st)];
                    left=left+length(f_st);
                end
                if u(4) == 1
                    Rep=[Rep; left+1 left+length(f_st)];
                    left=left+length(f_st);
                end
                if pos(j)+1>length(f_cur)
                    f_cur=[];
                else
                    f_cur=f_cur(pos(j)+1:end);
                end
                break;
            end                        
        end
    end
end

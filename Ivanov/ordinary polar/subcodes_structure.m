function [R0,R1,SPC,Rep] = subcodes_structure(f)
    N=length(f);
    f_cur=f;
    l_past=N;
    left=0;
    R0=[];
    R1=[];
    SPC=[];
    Rep=[];
    while l_past>0
        e=1;
        f_st=f_cur(1);
        while e>0
            while sum(ispattern(f_st))==1 && length(f_st)<length(f_cur)
                %if length(f_cur)>=2*length(f_st)
                    f_st=f_cur(1:2*length(f_st));
                %end                
            end            
            if sum(ispattern(f_st))==1 && length(f_st)>=length(f_cur)
                e=0;
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
            end
            if sum(ispattern(f_st))==0 && length(f_st)<=length(f_cur)
                while sum(ispattern(f_st))==0
                    f_st=f_st(1:length(f_st)/2);
                end
                u=ispattern(f_st);
                e=0;
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
            end
        end
        f_cur=f_cur(left+1:end);
        l_past=length(f_cur);
    end
end



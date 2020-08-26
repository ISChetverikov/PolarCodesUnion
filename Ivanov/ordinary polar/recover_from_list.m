function code_words=recover_from_list(list_u,prev)
if size(list_u,2)>1
    words_0=recover_from_list(list_u(:,1:end-1),prev(:,1:end-1));
    index=prev(:,end)>0;
    code_words=[words_0(prev(index,end),:) list_u(index,end)];
else
    code_words=list_u;
end
end
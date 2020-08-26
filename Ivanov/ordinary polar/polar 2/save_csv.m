function save_csv(fname)
    import_ws EsNos Errs Early
    f = fopen(fname, 'w');
    fprintf(f, 'esno,err,fail,frames,early\n');
    fclose(f);
    ii = (Errs(:,4) + Errs(:,3)) > 10;
    dlmwrite(fname, [EsNos(ii)', Errs(ii,3:5), Early(ii)], '-append');
end

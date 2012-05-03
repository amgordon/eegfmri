for s = 5:8%length(sa_all)
    
    
    par = EF_Params(sa_all{s});
    if par.goodSub
%         EF_wholeshebang(par, 'i');
%         par = EF_Params(sa_new{s});
%         EF_wholeshebang(par, 's');
%         par = EF_Params(sa_new{s});
        EF_wholeshebang(par, 'lcn');
        par = EF_Params(sa_all{s});
        EF_wholeshebang(par, 'h');
        EF_wholeshebang(par, 'z');
    end
    %hasHires = exist(par.hiresimg);
    %
    %EF_EEG_ResponseFunction(sa_good{s});
    %par = EF_Params(sa{s});
    %EF_wholeshebang(par, 'h');
    %EF_makeAnas(par)
end



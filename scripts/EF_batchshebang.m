for s = 2:length(sa_new)
    
    
    par = EF_Params(sa_new{s});
    if par.goodSub
%         EF_wholeshebang(par, 'i');
%         par = EF_Params(sa_new{s});
%         EF_wholeshebang(par, 's');
%         par = EF_Params(sa_new{s});
        EF_wholeshebang(par, 'rpe');
    end
    %hasHires = exist(par.hiresimg);
    %
    %EF_EEG_ResponseFunction(sa_good{s});
    %par = EF_Params(sa{s});
    %EF_wholeshebang(par, 'h');
    %EF_makeAnas(par)
end



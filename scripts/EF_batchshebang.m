for s = 2:8%length(sa_good)
    
    
    par = EF_Params(sa_good{s});
    if par.goodSub
        EF_wholeshebang(par, 'rpe');
    end
    %hasHires = exist(par.hiresimg);
    %
    %EF_EEG_ResponseFunction(sa_good{s});
    %par = EF_Params(sa{s});
    %EF_wholeshebang(par, 'h');
    %EF_makeAnas(par)
end






for i = 1:length(sa_good2)
    par = EF_Params(sa_good2{i});
    
    cd(par.analysisdir);
    
    
    delete('spmT*')
    delete('con*')
    
    load('SPM');
    SPM.xCon = [];
    
    save SPM SPM
    
end




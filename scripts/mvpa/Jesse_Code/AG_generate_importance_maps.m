function AG_generate_importance_maps(subj, results, results_IW, S)

subj = AG_interpret_weights_NEW(subj, results,results_IW, S);

vol_info = S.vol_info.vol_info;
voxel_inds = find(subj.masks{end}.mat);

for p = 1:length(subj.patterns)
    patNames{p} = subj.patterns{p}.name;
end

for it = 1:length(S.impType)
    for cds = 1:length(S.condnames)
        
        thisPat = strmatch([ 'impmap_' S.impType{it} '_' S.condnames{cds}] , patNames);
        
        impmap.(S.impType{it}).(S.condnames{cds}) = zeros(vol_info.dim);
        
        for rns = length(thisPat)
            
            
            
            
            
            
            temp{it}.cnd{cds} = zeros(vol_info.dim);
            
            temp{it}.cnd{cds}(voxel_inds) = subj.patterns{thisPat(rns)}.mat;
            
            impmap.(S.impType{it}).(S.condnames{cds}) = impmap.(S.impType{it}).(S.condnames{cds}) + temp{it}.cnd{cds};
            
        end  
        
            impmap.(S.impType{it}).(S.condnames{cds}) = impmap.(S.impType{it}).(S.condnames{cds})/length(thisPat);
            
            vol_info.dir = [S.importance_maps_dir '/' S.impType{it}];
            vol_info.fname = [ vol_info.dir '/' S.subj_id '_' S.impType{it} '_' S.condnames{cds} '.img'];
            
            if isempty(dir([vol_info.dir]))
                mkdir(vol_info.dir);
            end
            
            spm_write_vol(vol_info,impmap.(S.impType{it}).(S.condnames{cds}));
            
        
    end
end




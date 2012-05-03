function EF_ArtScansToOns(subpar)

% ---load par params if need be---
if isstruct(subpar) % if it is par_params struct
    par = subpar;
else % assume subject string
    par = AG1Params(subpar);
end


[zScoreA_cell, delta_cell] = EF_ArtRep_find_artifact_timepoints(subpar);

art.raw = load(fullfile(par.artrepdir, ['art_global_modified_' par.substr]));

for i=1:length(art.raw.zscoreA_cell)
    art.mot{i} = art.raw.delta_cell{i};
    art.sig{i} = art.raw.zscoreA_cell{i};
    
    idx.art.mot{i} = (art.mot{i} < par.art.motThresh);
    idx.art.sig{i} = (art.sig{i} <  par.art.sigThresh);
    
    idx.art.allNoArt{i} = idx.art.mot{i} .* idx.art.sig{i};
end

cd (par.artrepdir)
save ( 'ArtIDX', 'idx');






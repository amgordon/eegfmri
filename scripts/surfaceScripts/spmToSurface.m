function [ output_args ] = spmVolToSurface(sa, flags)
%given a subject with a full first-level analysis in volume space, project
%the relevant contrast maps onto a group level surface.




%% single subject analyses
for i=1:length(sa)
    
    par = EF_Params(sa);
    
    mapToPutOnSurface_h = dir(fullfile(par.analysisdir, 'con*.img'));
    for i=1:length(mapToPutOnSurface_h)
        par.mapToPutOnSurface{i} = fullfile(par.analdir, mapToPutOnSurface_h(i).name);
        
        par.registeredSurface{i} = fullfile(par.analdir, ['surf_' mapToPutOnSurface_h(i).name]);
        par.registeredSurfaceSmoothed{i} = fullfile(par.analdir, ['ssurf_' mapToPutOnSurface_h(i).name]);
    end
    
    %%
    
    % create surface from an anatomical volume
    reconCommand = ['recon-all -i ' par.hiresimg ' -s ' par.substr ' -all'];
    if ismember('r', flags), unix (reconCommand); end
    
    % register a functional volume to the 256 x 256 x 256 1mm anatomical images.
    % Output register.dat parameter file
    
    
    for j=1:length(par.mapToPutOnSurface)

        registerCommand = ['bbregister --s ' par.substr ' --mov ' par.mapToPutOnSurface ' --init-fsl --reg ' par.registerFile ' --bold'];
        if ismember('b', flags), unix (registerCommand); end
        
        for h=1:length(par.hemis)
            %put a con volume onto a group surface
            vol2SurfCommand = ['mri_vol2surf --mov ' par.mapToPutOnSurface ' --reg ' par.registerFile ' --hemi ' par.hemis{h} ' --trgsubject fsaverage --projfrac avg --o ' par.registeredSurface];
            if ismember('v', flags), unix (vol2SurfCommand); end
            
            % smooth statistical map
            smoothCommand = ['mri_surf2surf --sval ' par.registeredSurface ' --s ' par.substr ' --fwhm ' par.fs_smooth ' --tval ' par.registeredSurfaceSmoothed ' --hemi ' par.hemis{h}];
            if ismember('s', flags), unix (smoothCommand); end
            
            gpar.contrast(j).hemi(h).sub{s} = [par.registeredSurfaceSmoothed ' '];
        end
    end
end

%% group level analyses
%cd /Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/ef_091511/functional/scan01/ascan01_0006.nii
thisSPM = load(par.analysisDir, 'SPM.mat');
for j=1:length(gpar.contrast)
    for h=1:length(gpar.contrast(j).hemi)
        
        subsConcat = horzcat(gpar.contrast(j).hemi(h).sub{:});
        groupMapName = thisSPM.xCon(j).name;
        groupMapFile = fullfile(par.fs_anatdir, groupMapName);
        concatCommand = ['mri_concat ' subsConcat ' --o ' fullfile(groupMapFile, [par.hemis{h} '.' groupMapName]) ];
        unix(concatCommand);
        
        mri_glmfit --y lh.allsubs.avg.mgz --osgm --glmdir . --surf fsaverage lh --cortex
        
    end
end

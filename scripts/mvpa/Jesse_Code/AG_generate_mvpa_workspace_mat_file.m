function [subj] = AG_generate_mvpa_workspace_mat_file(subj_id, roi_name, img_files, mvpa_dir, runs_vector, exp_name, roi_file)

%exp_name = 'PAST';
%roi_file = ['/Users/alangordon/Accumulator_fMRI/AG1/fmri_data/ag1_021509/mvpa/onlyPPA.hdr'];

% load user-created filename and onsets lists into workspace

%load([mvpa_dir '/' data_imgs_to_use]); %loads predefined cell array called raw_filenames into memory

%num_runs = length(runs_vector); %calculate the number of runs (allows script to work flexibly for subjects with missing runs)
[subj] = AG_mvpa_load_and_preprocess_raw_data(subj_id, exp_name, roi_name, roi_file, img_files,  runs_vector);
% voxel_means = mean(subj.patterns{end}.mat,2);
% zero_voxel_inds = find(voxel_means==0);
% subj.patterns{end}.mat(zero_voxel_inds,:)=[]; % remove voxels with no data
% mask_vox = find(subj.masks{end}.mat);
% subj.masks{end}.mat(mask_vox(zero_voxel_inds))=0; % remove these voxels from mask as well
% display([num2str(length(zero_voxel_inds)) ' zero-value voxels removed from pattern and mask structs'])
% clear voxel_means zero_voxel_inds;
save_cmd = ['save ' mvpa_dir '/' subj_id '_' roi_name '.mat'];
eval(save_cmd);

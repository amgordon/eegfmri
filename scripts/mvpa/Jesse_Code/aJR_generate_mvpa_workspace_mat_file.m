function [subj num_runs num_TP_per_run] = aJR_generate_mvpa_workspace_mat_file(subj_id, roi_name, data_imgs_to_use, mvpa_dir)

par = par_params(subj_id);

exp_name = 'RRI_scan';
roi_file = ['/Users/bkuhl/research_stuff/experiments/fMRI/RRI_scan/MVPAarpeet/' roi_name '.img'];
% roi_file = ['C:/MVPA/mvpa/' roi_name '.img'];

num_TP_per_run = par.numvols; % AS: replaced with a vector of TP per run

% load user-created filename and onsets lists into workspace

% load([mvpa_dir '/' data_imgs_to_use]); %loads predefined cell array called raw_filenames into memory
num_runs = length(num_TP_per_run); %calculate the number of runs (allows script to work flexibly for subjects with missing runs)
                        % AS: replaced with the length of vectors
[subj] = aJR_mvpa_load_and_preprocess_raw_data(subj_id, exp_name, roi_name, roi_file, data_imgs_to_use, num_runs, num_TP_per_run);
% voxel_means = mean(subj.patterns{end}.mat,2);
% zero_voxel_inds = find(voxel_means==0);
% subj.patterns{end}.mat(zero_voxel_inds,:)=[]; % remove voxels with no data
% mask_vox = find(subj.masks{end}.mat);
% subj.masks{end}.mat(mask_vox(zero_voxel_inds))=0; % remove these voxels from mask as well
% display([num2str(length(zero_voxel_inds)) ' zero-value voxels removed from pattern and mask structs'])
% clear voxel_means zero_voxel_inds;
save_cmd = ['save ' mvpa_dir '/mvpa_results/subjworkspace.mat'];
eval(save_cmd);

function [marsS] = EF_marsbar_extract_betas()

roi_dir = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/group_analyses/analysis_group_by_RHitsVsCRs_LPI_400To700/RecVsCRsHC';

roi_names = {'AnGMaskedSphere1_-37_-69_30_roi.mat' 'AnGMaskedSphere2_-39_-74_47_roi.mat' 'AnGMaskedSphere3_-48_-52_24_roi.mat'...
'IPSMaskedSphere1_-34_-51_39_roi.mat' 'IPSMaskedSphere2_-36_-60_42_roi.mat' 'IPSMaskedSphere3_-36_-54_48_roi.mat'};

spm_name{1} = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/group_analyses/analysis_group_by_RHitsVsCRs_LPS_400To700/RecVsCRsHC/SPM.mat';
spm_name{2} = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/group_analyses/analysis_group_by_RHitsVsCRs_LPI_400To700/RecVsCRsHC/SPM.mat';
spm_name{3} = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/group_analyses/analysis_group_by_RHitsVsCRs_LPS_500To800/RecVsCRsHC/SPM.mat';
spm_name{4} = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/group_analyses/analysis_group_by_RHitsVsCRs_LPI_500To800/RecVsCRsHC/SPM.mat';
spm_name{5} = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/group_analyses/analysis_byIndDiff_LFS0p3to0p5sec/RecVsCRsHC/SPM.mat';
spm_name{6} = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/group_analyses/analysis_byIndDiff_LPS0p4to0p8sec/RecVsCRsHC/SPM.mat';
spm_name{7} = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/group_analyses/analysis_byIndDiff_MedPar0p4to0p8sec/RecVsCRsHC/SPM.mat';
spm_name{8} = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/group_analyses/analysis_byIndDiff_MedPar0p4to0p8sec/hits_HC_vs_LC/SPM.mat';
%spm_name{1} = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/group_analyses/analysis_hitsAndCRsByConf_Rec/RecVsCRsHC/SPM.mat';

for i=1:length(spm_name)
    
    thisSPM = load(spm_name{i});
    
    D = mardo(spm_name{i});
    %D = autocorr(D,'fmristat',2);
    
    for curr_clust=1:length(roi_names) % go through the list of clusters
        
        clusters_h =  [roi_dir '/' roi_names{curr_clust}];
        
        % Make marsbar design object
        roi_file = clusters_h;
        
        % Make marsbar ROI object
        R  = maroi(roi_file);
        
        % Fetch data into marsbar data object
        Y  = get_marsy(R, D, 'mean');
        
        % Get contrasts from original design
        xCon = get_contrasts(D);
        
        % Estimate design on ROI data
        E = estimate(D, Y);
        
        % Put contrasts from original design back into design object
        E = set_contrasts(E, xCon);
        
        % get design betas
        b = betas(E);
        
        marsS.clust(curr_clust).raw(i,:) = summary_data(Y);
        
        % get stats and stuff for all contrasts into statistics structure
        marsS.clust(curr_clust).spm(i) = compute_contrasts(E, 1:length(xCon));
    end
end



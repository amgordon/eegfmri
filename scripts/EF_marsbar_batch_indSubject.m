function [subjects_tc IPSC thisCon] = EF_marsbar_batch_indSubject(subjects)


conds = { 'Rem_Cor' 
    'CRs_HC'  
    };

roi_dir = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/group_analyses/analysis_group_by_RHitsVsCRs_LPI_400To700/RecVsCRsHC';
roi_names = {'AnG_sphere_roi.mat'	'IPS_sphere_roi.mat' 'AnG_sphere_roi2_roi.mat'	'IPS_sphere_roi2_roi.mat'  'AnG_sphere_roi3_roi.mat'	'IPS_sphere_roi3_roi.mat'};

thisCon = zeros(length(roi_names), length(subjects));
for curr_subj=1:length(subjects), % go through the list of subjects
    fprintf(1,'extracting TC from subject %s\n',subjects{curr_subj});
    spm_name = ['/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/' subjects{curr_subj} '/analysis_hitsAndCRsByConf_Rec/SPM.mat'];
    %spm_name = ['/Users/alangordon/Studies/AG1/Accumulator_fMRI/AG1/fmri_data2/' subjects{curr_subj} '/AD_balancedTrainingSet_50RandomIterations/SPM.mat'];
    D = mardo(spm_name);
    D = autocorr(D,'fmristat',2);
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
        % E = set_contrasts(E, xCon);
        
        % get design betas
        b = betas(E);
        
        thisCon(curr_clust, curr_subj) = b(5) - b(2);
        % get stats and stuff for all contrasts into statistics structure
        %marsS = compute_contrasts(E, 1:length(xCon));
        
        % Get definitions of all events in model
        [e_specs, e_names] = event_specs(E);
        n_events = size(e_specs, 2);
        
        % Bin size in seconds for FIR
        % bin_size = tr(E); %
        bin_size = 2;
        
        % Length of FIR in seconds
        fir_length = 20;
        
        % Number of FIR time bins to cover length of FIR
        bin_no = fir_length / bin_size;
        
        % Options - here 'single' FIR model, return estimated % signal change
        opts = struct('single', 1, 'percent', 1);
        
        % Return time courses for all events in fir_tc matrix
        
        %      for e_s = 1:n_events
        %        fir_tc(:, e_s) = event_fitted_fir(E, e_specs(:,e_s), bin_size, ...
        %                                          bin_no, opts);
        %      end
        
        %for cds = 1:length(condType)
            %subjects_tc(curr_subj).rois(curr_clust).conds = condType{cds};
            %subjects_tc(curr_subj).rois(curr_clust).condNames = conds{cds};
            clear fir_tc;
            for c = 1:length(conds)
                e_s = find(strcmp(conds{c}, e_names));
                
                fir_tc(:, c) = event_fitted_fir(E, e_specs(:,e_s), bin_size, bin_no, opts);
                
                [maxC, maxI] = max(fir_tc(2:((fir_length/bin_size)-1),c));
                
                subjects_tc(curr_subj).rois(curr_clust).conds = conds;
                
                subjects_tc(curr_subj).rois(curr_clust).roiName = roi_names{curr_clust};
                %subjects_tc(curr_subj).rois(curr_clust).maxval(c) = maxC;
                
                subjects_tc(curr_subj).rois(curr_clust).maxval(c) = maxC;
                neighbors = intersect(1:size(fir_tc,2), [maxI-1: maxI+1]);
                subjects_tc(curr_subj).rois(curr_clust).maxWithNeighbors(c) = mean(fir_tc(neighbors,c));
                
                subjects_tc(curr_subj).rois(curr_clust).raw = fir_tc;
            end
        %end
        
        
    end % end clusters loop
    
    subjects_tc(curr_subj).conds = conds;
    subjects_tc(curr_subj).e_specs = e_specs;
    subjects_tc(curr_subj).e_names = e_names;
    subjects_tc(curr_subj).n_events = n_events;
    
    
end % end subjects loop


%%
%interpreter information
TRs = 3:4;
conds = subjects_tc(1).rois(1).conds;

clear w IPSC mean_IPSC se_IPSC;
for i = [1:27]%length(subjects_tc)
    for j = 1:length(subjects_tc(i).rois)
    w{i,j}= subjects_tc(i).rois(j).raw;
    
    end
    
end

for k =1:size(w,2)
    w2 = cat(3, w{:,k});
    mean_tc_across_subs{k} = mean(w2, 3);
    
    IPSC{k} = squeeze(sum(w2(TRs,:,:),1))';
    
    mean_IPSC{k} = mean(IPSC{k});
    
    se_S_IPSC{k} = std(IPSC{k}) / sqrt(18);
    
    x1 = IPSC{k}(:);
    x2 = repmat([1:(length(subjects_tc)-1)]', length(conds),1);%repmat([1:length(subjects_tc)]', length(conds),1);
    x3 = sort(repmat([1:length(conds)]', (length(subjects_tc)-1), 1)); %sort(repmat([1:length(conds)]', length(subjects_tc), 1));
    
    X = [x1 x2 x3];
    
    aovRes{k} = ar_rmanova1(X);
    
    se_SXC_IPSC(k) = sqrt(aovRes{k}.MSE) / sqrt(18); %the more appropriate SE term to use!
    %figure;
   
    %%extra stuff pertaining to within and across group differences: clean it up
%     [p, STATS] = vartestn(IPSC{k}(:,[1 3]));
%     
%     ttest(diff(IPSC{k}(:,1:2)'), diff(IPSC{k}(:,2:3)'))
    
    
    %plot(mean_tc_across_subs{k});
end
 a = cat(3,IPSC{:});
 b = squeeze(mean(a,1));
c = vertcat(se_S_IPSC{:})';

for i = [1 5]
 figure; barweb(b(i:(i+3),:), c(i:(i+3),:)' );
end                         
                                   

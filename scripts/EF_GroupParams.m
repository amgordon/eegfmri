function gpar = EF_GroupParams(subjArray)
gpar.subjArray = subjArray;

gpar.conGroup = { 'analysis_hitsVsCRs_by_LPIAmp' };

gpar.tasks = {'analysis_hitsVsCRs_by_LPIAmp_12subs'};

%gpar.expt_dir = '/Users/alangordon/mounts/w5/alan/perceptMnemonic/fmri_data/';
gpar.expt_dir = '/biac4/wagner/biac3/wagner5/alan/eegfmri/fmri_data/';

gpar.modelTemplate = '/Users/gordonam/Studies/AG1/scripts/GroupTemplate.mat';

gpar.constat = 'T';
gpar.exMask = [];


for t = 1:length(gpar.tasks)
    
    gpar.task{t}.conTemplate =     fullfile(gpar.expt_dir, 'ef_042912', gpar.conGroup{t}, 'SPM');
    
    ldTemp = load(gpar.task{t}.conTemplate);
    
    gpar.task{t}.SPMcons = ldTemp.SPM.xCon;
    
    gpar.nCovs = 0;
    
    %erp_cov = load('/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/erp_data/PMMeanZscore_0p4to0p8sec_RHCH_LCH_LCCR_HCCR.mat');
 
    %for i=1:length(subjArray), subIDs(i) = EF_num2Sub(subjArray{i}); end
    %idxSubsToInclude = ismember(erp_cov.Y(:,1), subIDs);
    
    %gpar.covVec = erp_cov.Y(idxSubsToInclude,3) - erp_cov.Y(idxSubsToInclude,2);
    %gpar.covName = 'subjectWise_LCHits_vs_HCHitsAndRHits_MedParietal_400To800';
    gpar.covVec = [];
    gpar.covName = [];
     
    for c= 1:length(gpar.task{t}.SPMcons);
        gpar.exMask = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/groupMask/inclusive_mask.img';
        
        gpar.task{t}.cons{c}.dir = {fullfile(gpar.expt_dir, 'group_analyses', gpar.tasks{t},gpar.task{t}.SPMcons(c).name)};
        gpar.task{t}.cons{c}.name = gpar.task{t}.SPMcons(c).name;
        
        for s = 1:length(gpar.subjArray)
            gpar.task{t}.cons{c}.scans{s} = fullfile(gpar.expt_dir, gpar.subjArray{s}, gpar.conGroup{t}, ['con_' prepend(num2str(c), 4) '.img']);
        end
        
    end
end
function [] = EF_EEG_EstimateIndividualBetas(par)

orig_par = par;
[~, idx] = EF_BehAnalyzer(par);

for i=1:length(idx.allTrials) 
    
    thisSess = idx.sess(i);
    sessBaselineScans = sum(par.numvols(1:(thisSess-1)));
    sessBaselineTrials = sum(idx.sess<thisSess);
    allTrialsThisSess = idx.allTrials(idx.sess==thisSess) - sessBaselineScans;

    idxScansThisSess = 1+sum(par.numvols(1:(thisSess-1))):sum(par.numvols(1:thisSess));
    idxThisTrialWithinSess = i - sessBaselineTrials;
    
    par.swascanfiles = orig_par.swascanfiles(idxScansThisSess,:);
    par.numscans = orig_par.numscans(thisSess);
    
    onsets{1} = allTrialsThisSess(idxThisTrialWithinSess);
    onsets{2} = setdiff(allTrialsThisSess, allTrialsThisSess(idxThisTrialWithinSess));
    
    thisTrialName = ['trial' prepend(num2str(i),3)];
    names = {thisTrialName 'otherTrials'};
    
    durations = num2cell(zeros(size(names)));
   
    if ~exist(par.analysisdir,'dir')
        mkdir(par.analysisdir);
    end
    
    cd(par.analysisdir);
    save ons onsets names durations
    
    regs = [];
    save regs;

    
    EF_mod_spec(par);
    EF_mod_est(par);
    
    movefile('beta_0001.hdr', [thisTrialName '.hdr']);
    movefile('beta_0001.img', [thisTrialName '.img']);
    
    delete('beta*', 'ResMS*', 'SPM*', 'RPV*', 'mask*', 'ons*', 'regs*');
end


end


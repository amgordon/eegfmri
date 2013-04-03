function [subj results] = AG_extractMeanSignal(subj, S)

%subj = load_spm_mask(subj,S.ROI1_name,S.ROI1_file);
%subj = load_spm_mask(subj,S.ROI2_name,S.ROI2_file);

selectorMat = squeeze(get_group_as_matrix(subj,'selector',S.thisSigIntenseSelector));

if S.xval
    train_idx = (selectorMat~=0);
    test_idx = (selectorMat~=0);
else
    train_idx = (selectorMat==1);
    test_idx  = (selectorMat==2);
end

%subj = create_pattern_from_mask(subj,S.preprocPatCondensedName,S.ROI1_name,S.ROI1PatName);
%subj = create_pattern_from_mask(subj,S.preprocPatCondensedName,S.ROI2_name,S.ROI2PatName);

if S.defineROIsFromANOVAFS
    
    if ~isempty(strfind(S.thisSigIntenseSelector, 'xval'))
        sel = squeeze(get_group_as_matrix(subj,'selector',S.thisSigIntenseSelector));
    else
        sel = get_mat(subj,'selector',S.thisSigIntenseSelector);
    end
    
    sROI1 = get_mat(subj,'pattern',S.preprocPatCondensedName);
    sROI2 = get_mat(subj,'pattern',S.preprocPatCondensedName);
    
    reg = (get_mat(subj,'regressors', 'conds'));
    
    if S.xval
        for j = 1:size(sel,1)
            
            thisTestIt = find(sel(j,:)==2);
            
            trainSig1 = sROI1(:,find((sel(j,:)==1) .* (reg(1,:)==1)));
            trainSig2 = sROI1(:,find((sel(j,:)==1) .* (reg(2,:)==1)));
            
            [a1 a2 a3 a4] = ttest2(trainSig1', trainSig2');
            sorted_tstats = sort(a4.tstat);
            
            thresh1 =  sorted_tstats(500);
            thresh2 =  sorted_tstats(end - 500);
            
            idxROI2 = find(a4.tstat<=thresh1);
            idxROI1 = find(a4.tstat>thresh2);
            
            testSig1(thisTestIt) = mean(sROI1(idxROI1,thisTestIt));
            testSig2(thisTestIt) = mean(sROI1(idxROI2,thisTestIt));
            
        end
        
        test_targs = reg;
    else
        thisTestIt = find(sel==2);
        
        trainSig1 = sROI1(:,find((sel==1) .* (reg(1,:)==1)));
        trainSig2 = sROI1(:,find((sel==1) .* (reg(2,:)==1)));
        
        [a1 a2 a3 a4] = ttest2(trainSig1', trainSig2');
        sorted_tstats = sort(a4.tstat);
        
        thresh1 =  sorted_tstats(500);
        thresh2 =  sorted_tstats(end - 500);
        
        idxROI2 = find(a4.tstat<=thresh1);
        idxROI1 = find(a4.tstat>thresh2);
        
        testSig1 = mean(sROI1(idxROI1,thisTestIt));
        testSig2 = mean(sROI1(idxROI2,thisTestIt));
        
        test_targs = reg(:,thisTestIt);
    end
    
else
    %sigROI1 = get_mat(subj,'pattern',S.ROI1PatName);
    %sigROI2 = get_mat(subj,'pattern',S.ROI2PatName);
    
    sigROI1 = get_mat(subj,'pattern',S.preprocPatCondensedName);
    sigROI2 = sigROI1;
        
    all_targs = get_mat(subj,'regressors','conds');
    
    train_targs = all_targs(:,train_idx);
    test_targs = all_targs(:,test_idx);
    
    
    DV = (2-train_targs(1,:)');
    IVs = [mean(sigROI1(:,train_idx))'  mean(sigROI2(:,train_idx))'];
    
    TestVals = [mean(sigROI1(:,test_idx))'  mean(sigROI2(:,test_idx))'];
    
    testSig1 = mean(sigROI1(:,test_idx));
    testSig2 = mean(sigROI2(:,test_idx));
end

if S.logreg_2Features %%%NOTE THAT THIS IS NOT YET IMPLEMENTED FOR N-FOLD XVALIDATION!!!   ONLY FOR TRAIN AND TEST ON DIFFERENT SETS!
    %run a logistic regression, using each of the mean signals as IVs
    Betas = mnrfit(IVs, DV);
    acts(1,:) = Betas(1) + TestVals * Betas(2:3);
    acts(2,:) = -acts(1,:);
elseif S.zscoreIntensityVals
    % just do comparisons of the mean signal value
    acts(1,:) = zscore(testSig1);
    acts(2,:) = zscore(testSig2);  
else
    acts(1,:) = testSig1;
    acts(2,:) = testSig2;
end

results.iterations(1).acts = acts;
perfmet = perfmet_maxclass(acts,test_targs,[]);

results.iterations(1).perfmet = perfmet;

results.total_perf = perfmet.perf;

fprintf( '\ngot total_perfs - %s \n', num2str(perfmet.perf));



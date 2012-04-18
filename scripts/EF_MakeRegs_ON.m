function EF_MakeRegs_ON(par)

[res idx] = EF_BehAnalyzer(par);

erpDat = load(par.erpFile);

%combine all trials into one 400 length vector??

% 
for i = 1:length(par.usedVols)
    erpidx.all_h{i} = erpDat.subj(par.subNo).run(i).all.goodtr;
    binsize = size(erpDat.subj(par.subNo).run(i).all.bins);
    
    erpidx.sigAll_h{i} = nan(binsize(1), binsize(2), size(erpidx.all_h{i},1));
    erpidx.sigAll_h{i}(:,:,erpidx.all_h{i}==1) = erpDat.subj(par.subNo).run(i).all.bins;
    
    erpidx.sigAll_mean{i} = squeeze(mean(mean(erpidx.sigAll_h{i}(:,6:8,:))));
end
% 
 erpidx.sigAll = vertcat(erpidx.sigAll_mean{:});
 erpidx.all = vertcat(erpidx.all_h{:})';
% 
if strcmp(par.substr, 'ef_072111')
    erpidx.all = [erpidx.all 0];
    erpidx.sigAll = [erpidx.sigAll; 0];
end

switch par.ONAnalysisType
    
    case 'hitsVsCRs'
        onsets{1} = idx.allTrials(find(idx.old .* idx.respOld));
        onsets{2} = idx.allTrials(find(idx.new .* idx.respNew));
        onsets{3} = idx.allTrials(find(idx.old .* ~idx.respOld + idx.new .* ~idx.respNew));
        names = {'hits' 'CRs' 'junk'};
    case 'hitsVsCRs_HC'
        onsets{1} = idx.allTrials(find(idx.old .* idx.respOld .* idx.highConf));
        onsets{2} = idx.allTrials(find(idx.new .* idx.respNew .* idx.highConf));
        onsets{3} = idx.allTrials(find(idx.old .* idx.respOld .* idx.lowConf));
        onsets{4} = idx.allTrials(find(idx.new .* idx.respNew .* idx.lowConf));
        onsets{5} = idx.allTrials(find(idx.old .* ~idx.respOld + idx.new .* ~idx.respNew));
        names = {'hits_HC' 'CRs_HC' 'hits_LC' 'CRs_LC' 'junk'};
    case 'hitsVsCRs_byERP'
        
        onsets{1} = idx.allTrials(find(erpidx.all .* idx.old .* idx.respOld ));
        onsets{2} = idx.allTrials(find(erpidx.all .* idx.new .* idx.respNew ));
        onsets{3} = idx.allTrials(find(erpidx.all .* idx.old .* idx.respNew ));
        onsets{4} = idx.allTrials(find(erpidx.all .* idx.new .* idx.respOld ));
        onsets{5} = setdiff(idx.allTrials(find(idx.old + idx.new)), [onsets{1:4}]);
        
        ERPAmp_all = erpidx.sigAll(find(erpidx.all .* idx.cor ));
        ERPAmp_Hits = erpidx.sigAll(find(erpidx.all .* idx.cor .* idx.old ));
        ERPAmp_CRs = erpidx.sigAll(find(erpidx.all .* idx.cor .* idx.new ));
        
        pmod(1).param = {ERPAmp_Hits};
        pmod(1).name = {'ERPAmp_Hits'};
        pmod(1).poly = {1};
        
        pmod(2).param = {ERPAmp_CRs};
        pmod(2).name = {'ERPAmp_CRs'};
        pmod(2).poly = {1};
        
        names = {'hits' 'CRs' 'misses' 'FAs' 'junk'};
        

    case 'hitsVsCRs_byERP_HC'
          
        onsets{1} = idx.allTrials(find(erpidx.all .* idx.old .* idx.respOld .* idx.highConf));
        onsets{2} = idx.allTrials(find(erpidx.all .* idx.new .* idx.respNew .* idx.highConf));
        onsets{3} = idx.allTrials(find(erpidx.all .* idx.old .* idx.respNew .* idx.highConf));
        onsets{4} = idx.allTrials(find(erpidx.all .* idx.new .* idx.respOld .* idx.highConf));
        onsets{5} = idx.allTrials(find(erpidx.all .* idx.old .* idx.respOld .* idx.lowConf));
        onsets{6} = idx.allTrials(find(erpidx.all .* idx.new .* idx.respNew .* idx.lowConf));
        onsets{7} = idx.allTrials(find(erpidx.all .* idx.old .* idx.respNew .* idx.lowConf));
        onsets{8} = idx.allTrials(find(erpidx.all .* idx.new .* idx.respOld .* idx.lowConf));
        onsets{9} = setdiff(idx.allTrials(find(idx.old + idx.new)), [onsets{1:4}]);
        
        ERPAmp_all = erpidx.sigAll(find(erpidx.all .* idx.cor ));
        ERPAmp_Hits = erpidx.sigAll(find(erpidx.all .* idx.cor .* idx.old ));
        ERPAmp_CRs = erpidx.sigAll(find(erpidx.all .* idx.cor .* idx.new ));
        
        pmod(1).param = {ERPAmp_Hits};
        pmod(1).name = {'ERPAmp_Hits'};
        pmod(1).poly = {1};
        
        pmod(2).param = {ERPAmp_CRs};
        pmod(2).name = {'ERPAmp_CRs'};
        pmod(2).poly = {1};
        
        names = {'hits_HC' 'CRs_HC' 'misses_HC' 'FAs_HC' 'hits_LC' 'CRs_LC' 'misses_LC' 'FAs_LC' 'junk'};
        
    case 'alltrials_byERP'
        
        
        onsets{1} = idx.allTrials(find(erpidx.all .* idx.cor ));
        onsets{2} = setdiff(idx.allTrials(idx.old+idx.new==1), [onsets{1}]);
        
        ERPAmp_all = erpidx.sigAll(find(erpidx.all .* idx.cor ));
        
        pmod(1).param = {ERPAmp_all};
        pmod(1).name = {'ERPAmp_all'};
        pmod(1).poly = {1};
        
        names = {'all_ERP' 'junk'};
        
    case 'alltrials_byERP_HC'

        
        onsets{1} = idx.allTrials(find(erpidx.all .* idx.cor .* idx.highConf));
        onsets{2} = setdiff(idx.allTrials(idx.old+idx.new==1), [onsets{1}]);
        
        ERPAmp_all = erpidx.sigAll(find(erpidx.all .* idx.cor .* idx.highConf));
        
        pmod(1).param = {ERPAmp_all};
        pmod(1).name = {'ERPAmp_all'};
        pmod(1).poly = {1};
        
        names = {'all_HC' 'junk'};
end

durations = num2cell(zeros(size(names)));

sessReg = zeros(sum(par.usedVols),length(par.usedVols)-1);
for i = 1:(length(par.usedVols)-1)
    sessReg(sum(par.usedVols(1:i-1))+1 : sum(par.usedVols(1:i)),i) = ones(par.usedVols(i),1);
end

R = sessReg;

if ~exist(par.analysisdir)
    mkdir(par.analysisdir);
end

cd (par.analysisdir);

onsets

%save ons onsets names durations pmod;
save ons onsets names durations pmod;
save regs R;




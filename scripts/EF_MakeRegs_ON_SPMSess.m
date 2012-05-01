function EF_MakeRegs_ON(par)

[res idx] = EF_BehAnalyzer(par);

%erpDat = load(par.erpFile);
%RTERPDat = load(par.rawERPFile);
%power = load(par.rawSpectralFile);

channel_rois;
%combine all trials into one 400 length vector??

% 
% for i = 1:length(par.usedVols)
%     erpidx.all_h{i} = erpDat.subj(par.subNo).run(i).all.goodtr;
%     binsize = size(erpDat.subj(par.subNo).run(i).all.bins);
%     
%     erpidx.sigAll_h{i} = nan(binsize(1), binsize(2), size(erpidx.all_h{i},1));
%     erpidx.sigAll_h{i}(:,:,erpidx.all_h{i}==1) = erpDat.subj(par.subNo).run(i).all.bins;
%     
%     erpidx.sigAll_mean{i} = squeeze(mean(mean(erpidx.sigAll_h{i}(:,6:8,:))));
% end
% % 
% %
% erpidx.sigAll = vertcat(erpidx.sigAll_mean{:});
% erpidx.all = vertcat(erpidx.all_h{:})';
% 

eeg_fmri_on_subject_info;
erpidx.all = ones(size(idx.allTrials));
erpidx.all(S.badtr{par.subNo}) = 0;


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
        
    case 'hitsVsCRs_ByEEGCertainty'
        
        if par.concatAcrossRuns
            idx.sessToSPM = ones(size(idx.sess));
        else
            idx.sessToSPM = idx.sess;
        end
        
        selectedEEGDat = load(fullfile(par.classmatDir, 'classifyLeftAnG.mat'));
        idxThisSub = strcmp(selectedEEGDat.subs, par.substr);
        thisS = selectedEEGDat.S{idxThisSub};
        eegCert = nan(size(thisS.idxCorMemory));
        eegCert(thisS.idxCorMemory) = selectedEEGDat.res{idxThisSub}.classifyMemoryStateWithEEG.mod.YProb;

        
        for i=1:length(par.sess)
            sess(i).onsets{1} = idx.onsets(find(thisS.idxCorMemory .* (idx.sessToSPM==i)'));
            sess(i).onsets{2} = idx.onsets(find(~thisS.idxCorMemory .* (idx.old .* idx.respOld + idx.new .* idx.respNew) .* (idx.sessToSPM==i)'));
            sess(i).onsets{3} = idx.onsets(find((idx.old .* ~idx.respOld + idx.new .* ~idx.respNew)  .* (idx.sessToSPM==i)'));
            
            sess(i).names = {'correctWithGoodERPs' 'correctWithBadERPs' 'incorrect'};
            
            % remove empty regressors
            toRemove = [];
            for j = 1:length(sess(i).onsets)
                if isempty(sess(i).onsets{j})
                    toRemove = [toRemove j];
                end
            end
            

            
            sess(i).pmod(1).param{1} = double(idx.old(thisS.idxCorMemory' & idx.sessToSPM==i)');
            sess(i).pmod(1).name{1} = 'hitsVsCRs';
            sess(i).pmod(1).poly{1} = 1;
            
            sess(i).pmod(1).param{2} = eegCert(thisS.idxCorMemory' & idx.sessToSPM==i)';
            sess(i).pmod(1).name{2} = 'correctByEEGCert';
            sess(i).pmod(1).poly{2} = 1;
            
            sess(i).onsets(toRemove) = [];
            sess(i).names(toRemove) = [];
    
        end
 
          

        

    case 'hitsVsCRs_HC'
        
        if par.concatAcrossRuns
            idx.sessToSPM = ones(size(idx.sess));
        else
            idx.sessToSPM = idx.sess;
        end
        
        for i=1:length(par.sess)
            sess(i).onsets{1} = idx.onsets(find(idx.old .* idx.respOld .* idx.highConf .* (idx.sessToSPM==i)'));
            sess(i).onsets{2} = idx.onsets(find(idx.new .* idx.respNew .* idx.highConf .* (idx.sessToSPM==i)'));
            sess(i).onsets{3} = idx.onsets(find(idx.old .* idx.respOld .* idx.lowConf .* (idx.sessToSPM==i)'));
            sess(i).onsets{4} = idx.onsets(find(idx.new .* idx.respNew .* idx.lowConf .* (idx.sessToSPM==i)'));
            sess(i).onsets{5} = idx.onsets(find((idx.old .* ~idx.respOld + idx.new .* ~idx.respNew)  .* (idx.sessToSPM==i)'));
            
            sess(i).names = {'hits_HC' 'CRs_HC' 'hits_LC' 'CRs_LC' 'junk'};
            
            % remove empty regressors
            toRemove = [];
            for j = 1:length(sess(i).onsets)
                if isempty(sess(i).onsets{j})
                    toRemove = [toRemove j];
                end
            end
            
            sess(i).onsets(toRemove) = [];
            sess(i).names(toRemove) = [];
        end
        
    case 'buttonPresses'
        idx.RTsExist = (RTERPDat.data.RTs~=2.1)';
        
        if strcmp(par.substr, 'ef_092711') % last 3 sessions have bad erp data.
            erpidx.all = erpidx.all(1:160);
            idx.RTsExist = idx.RTsExist(1:160);
        end
        
        onsets{1} = idx.allTrials(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'1') ));
        onsets{2} = idx.allTrials(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'2') ));
        onsets{3} = idx.allTrials(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'3') ));
        onsets{4} = idx.allTrials(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'4') ));
        onsets{5} = idx.allTrials(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'5') ));
        onsets{6} = idx.allTrials(find(idx.noResp .* ~idx.fix));
        onsets{7} = idx.allTrials(find((~idx.noResp .* ~idx.fix ) .* (~idx.RTsExist + ~erpidx.all)) );
        
        names = {'finger1' 'finger2' 'finger3' 'finger4' 'finger5' 'noResp' 'junk'};     
        
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
        
    case 'buttonPresses_byERP_RTLocked'
        
        %
        idx.RTsExist = (RTERPDat.data.RTs~=2.1)';
        
        if strcmp(par.substr, 'ef_092711') % last 3 sessions have bad erp data.
            erpidx.all = erpidx.all(1:160);
            idx.RTsExist = idx.RTsExist(1:160);
        end
        
        onsets{1} = idx.allTrials(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'1') ));
        onsets{2} = idx.allTrials(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'2') ));
        onsets{3} = idx.allTrials(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'3') ));
        onsets{4} = idx.allTrials(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'4') ));
        onsets{5} = idx.allTrials(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'5') ));
        onsets{6} = idx.allTrials(find(idx.noResp .* ~idx.fix));
        onsets{7} = idx.allTrials(find((~idx.noResp .* ~idx.fix ) .* (~idx.RTsExist + ~erpidx.all)) );
        
        %make sure erpidx.sigAll is RT locked
        erpAmp_h{1} = RTERPDat.data.trialdata(chnls.central.LCI, find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'1') ),par.criticalSpectralSamples);
        erpAmp_h{2} = RTERPDat.data.trialdata(chnls.central.LCI, find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'2') ),par.criticalSpectralSamples);
        erpAmp_h{3} = RTERPDat.data.trialdata(chnls.central.LCI, find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'3') ),par.criticalSpectralSamples);
        erpAmp_h{4} = RTERPDat.data.trialdata(chnls.central.LCI, find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'4') ),par.criticalSpectralSamples);
        erpAmp_h{5} = RTERPDat.data.trialdata(chnls.central.LCI, find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'5') ),par.criticalSpectralSamples);
        
        erpAmp{1} = mean(squeeze(mean(erpAmp_h{1})),2);
        erpAmp{2} = mean(squeeze(mean(erpAmp_h{2})),2);
        erpAmp{3} = mean(squeeze(mean(erpAmp_h{3})),2);
        erpAmp{4} = mean(squeeze(mean(erpAmp_h{4})),2);
        erpAmp{5} = mean(squeeze(mean(erpAmp_h{5})),2);
        
        pmod(1).param = {erpAmp{1}};
        pmod(1).name = {'ERPAmp_finger1'};
        pmod(1).poly = {1};
        
        pmod(2).param = {erpAmp{2}};
        pmod(2).name = {'ERPAmp_finger2'};
        pmod(2).poly = {1};
        
        pmod(3).param = {erpAmp{3}};
        pmod(3).name = {'ERPAmp_finger3'};
        pmod(3).poly = {1};
        
        pmod(4).param = {erpAmp{4}};
        pmod(4).name = {'ERPAmp_finger4'};
        pmod(4).poly = {1};
        
        pmod(5).param = {erpAmp{5}};
        pmod(5).name = {'ERPAmp_finger5'};
        pmod(5).poly = {1};
        
        pmod(6).param = [];
        pmod(6).name = [];
        pmod(6).poly = [];
        
        pmod(7).param = [];
        pmod(7).name = [];
        pmod(7).poly = [];
        
        names = {'finger1' 'finger2' 'finger3' 'finger4' 'finger5' 'noResp' 'junk'};
        
    case 'alltrials_byERP'
        
        idxCorTrials = find(erpidx.all .* idx.cor );
        
        onsets{1} = idx.allTrials(idxCorTrials);
        onsets{2} = setdiff(idx.allTrials(idx.old+idx.new==1), [onsets{1}]);
        
        selectedERPDat = RTERPDat.data.trialdata(chnls.central.LCI, idxCorTrials, par.criticalERPSamples);
        meanERPAcrossROI = squeeze(mean(selectedERPDat,1));
        meanERPAcrossROIAndTime = squeeze(mean(meanERPAcrossROI,2));
        
        ERPAmp_all = meanERPAcrossROIAndTime;
        
        pmod(1).param = {ERPAmp_all};
        pmod(1).name = {'ERPAmp_all'};
        pmod(1).poly = {1};
        
        names = {'all_ERP' 'junk'};
        
    case 'alltrials_byERP_HC'
        
        idxHCCorTrials = find(erpidx.all .* idx.cor .* idx.highConf);
        
        onsets{1} = idx.allTrials(idxHCCorTrials);
        onsets{2} = setdiff(idx.allTrials(idx.old+idx.new==1), [onsets{1}]);
        
        selectedERPDat = RTERPDat.data.trialdata(chnls.parietal.LPS, idxHCCorTrials, 350:450);
        meanERPAcrossROI = squeeze(mean(selectedERPDat,1));
        meanERPAcrossROIAndTime = squeeze(mean(meanERPAcrossROI,2));
        
        ERPAmp_all = meanERPAcrossROIAndTime;
        
        pmod(1).param = {ERPAmp_all};
        pmod(1).name = {'ERPAmp_all'};
        pmod(1).poly = {1};
        
        names = {'all_HC' 'junk'};
        
    case 'alltrials_by_spectralPower'
        
        idxCorTrials = find(erpidx.all .* idx.cor );
        
        onsets{1} = idx.allTrials(idxCorTrials);
        onsets{2} = setdiff(idx.allTrials(idx.old+idx.new==1), [onsets{1}]);
        
        selectedERPDat = power.data.normpower(chnls.central.LCI, par.bandOfInterest, idxCorTrials, par.criticalSpectralSamples);
        meanERPAcrossROI = squeeze(mean(selectedERPDat,1));
        meanERPAcrossROIAndTime = squeeze(mean(meanERPAcrossROI,2));
        
        ERPAmp_all = meanERPAcrossROIAndTime;
        
        pmod(1).param = {ERPAmp_all};
        pmod(1).name = {'ERPAmp_all'};
        pmod(1).poly = {1};
        
        names = {'all_ERP' 'junk'};
        
    case 'buttonPress_by_spectralPower'
                 
        onsets{1} = idx.allTrials(find(erpidx.all .* strcmp(idx.firstResp,'1') ));
        onsets{2} = idx.allTrials(find(erpidx.all .* strcmp(idx.firstResp,'2') ));
        onsets{3} = idx.allTrials(find(erpidx.all .* strcmp(idx.firstResp,'3') ));
        onsets{4} = idx.allTrials(find(erpidx.all .* strcmp(idx.firstResp,'4') ));
        onsets{5} = idx.allTrials(find(erpidx.all .* strcmp(idx.firstResp,'5') ));

        selectedERPDat = power.data.normpower(chnls.central.LCI, par.bandOfInterest, :, par.criticalSpectralSamples);
        meanERPAcrossROI = squeeze(mean(selectedERPDat,1));
        meanERPAcrossROIAndTime = squeeze(mean(meanERPAcrossROI,2));
        
        ERPAmp_all = meanERPAcrossROIAndTime;
        
        erpAmp{1} = ERPAmp_all(find(erpidx.all .* strcmp(idx.firstResp,'1') ));
        erpAmp{2} = ERPAmp_all(find(erpidx.all .* strcmp(idx.firstResp,'2') ));
        erpAmp{3} = ERPAmp_all(find(erpidx.all .* strcmp(idx.firstResp,'3') ));
        erpAmp{4} = ERPAmp_all(find(erpidx.all .* strcmp(idx.firstResp,'4') ));
        erpAmp{5} = ERPAmp_all(find(erpidx.all .* strcmp(idx.firstResp,'5') ));
        
        pmod(1).param = {erpAmp{1}};
        pmod(1).name = {'ERPAmp_finger1'};
        pmod(1).poly = {1};
        
        pmod(2).param = {erpAmp{2}};
        pmod(2).name = {'ERPAmp_finger2'};
        pmod(2).poly = {1};
        
        pmod(3).param = {erpAmp{3}};
        pmod(3).name = {'ERPAmp_finger3'};
        pmod(3).poly = {1};
        
        pmod(4).param = {erpAmp{4}};
        pmod(4).name = {'ERPAmp_finger4'};
        pmod(4).poly = {1};
        
        pmod(5).param = {erpAmp{5}};
        pmod(5).name = {'ERPAmp_finger5'};
        pmod(5).poly = {1};
        
        names = {'finger1' 'finger2' 'finger3' 'finger4' 'finger5'}; 
        
    case 'buttonPress_by_spectralPower_RTLocked'
        
        idx.RTsExist = (RTERPDat.data.RTs~=2.1)';
        
        if strcmp(par.substr, 'ef_092711') % last 3 sessions have bad erp data.
            erpidx.all = erpidx.all(1:160);
            idx.RTsExist = idx.RTsExist(1:160);
        end

        onsets{1} = idx.allTrials(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'1') ));
        onsets{2} = idx.allTrials(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'2') ));
        onsets{3} = idx.allTrials(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'3') ));
        onsets{4} = idx.allTrials(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'4') ));
        onsets{5} = idx.allTrials(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'5') ));
        onsets{6} = idx.allTrials(find(idx.noResp .* ~idx.fix));
        onsets{7} = idx.allTrials(find((~idx.noResp .* ~idx.fix ) .* (~idx.RTsExist + ~erpidx.all)) );

        
        selectedERPDat = power.data.normpower(chnls.central.LCI, par.bandOfInterest, :, par.criticalSpectralSamples);
        meanERPAcrossROI = squeeze(mean(selectedERPDat,1));
        meanERPAcrossROIAndTime = squeeze(mean(meanERPAcrossROI,2));
        
        ERPAmp_all = meanERPAcrossROIAndTime;
        
        erpAmp{1} = ERPAmp_all(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'1') ));
        erpAmp{2} = ERPAmp_all(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'2') ));
        erpAmp{3} = ERPAmp_all(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'3') ));
        erpAmp{4} = ERPAmp_all(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'4') ));
        erpAmp{5} = ERPAmp_all(find(erpidx.all .* idx.RTsExist .* strcmp(idx.firstResp,'5') ));
        
        pmod(1).param = {erpAmp{1}};
        pmod(1).name = {'ERPAmp_finger1'};
        pmod(1).poly = {1};
        
        pmod(2).param = {erpAmp{2}};
        pmod(2).name = {'ERPAmp_finger2'};
        pmod(2).poly = {1};
        
        pmod(3).param = {erpAmp{3}};
        pmod(3).name = {'ERPAmp_finger3'};
        pmod(3).poly = {1};
        
        pmod(4).param = {erpAmp{4}};
        pmod(4).name = {'ERPAmp_finger4'};
        pmod(4).poly = {1};
        
        pmod(5).param = {erpAmp{5}};
        pmod(5).name = {'ERPAmp_finger5'};
        pmod(5).poly = {1};
        
        pmod(6).param = [];
        pmod(6).name = [];
        pmod(6).poly = [];
        
        pmod(7).param = [];
        pmod(7).name = [];
        pmod(7).poly = [];
        
        names = {'finger1' 'finger2' 'finger3' 'finger4' 'finger5' 'noResp' 'junk'};
        
end

if ~exist(par.analysisdir)
    mkdir(par.analysisdir);
end

cd (par.analysisdir);

sessReg = zeros(sum(par.usedVols),length(par.usedVols)-1);
for i = 1:(length(par.usedVols)-1)
    sessReg(sum(par.usedVols(1:i-1))+1 : sum(par.usedVols(1:i)),i) = ones(par.usedVols(i),1);
end


% 
% if par.concatAcrossRuns
%     allOns = [onsets{:}];
%     
%     if (length(allOns)~=length(unique(allOns)))
%         error ('repeated onsets');
%     end
%     
%     durations = num2cell(zeros(size(names)));
%     
%     
%     %save ons onsets names durations pmod;
%     if exist('pmod')
%         save ons onsets names durations pmod;
%     else
%         save ons onsets names durations;
%     end
%     R = sessReg;
%     save regs R;
% else
%

for i=1:length(par.sess)
    
    onsets = sess(i).onsets;
    names = sess(i).names;
    durations = num2cell(zeros(size(names)));
    pmod = sess(i).pmod;
    
    R = [];
    
    if exist('pmod')
        save(sprintf('ons%g',i), 'onsets', 'names', 'durations', 'pmod');
    else
        save(sprintf('ons%g',i), 'onsets', 'names', 'durations');
    end
    
    save(sprintf('regs%g', i), 'R');
end
% end






      


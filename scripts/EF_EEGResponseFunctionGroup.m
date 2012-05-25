function [groupDat res S] = EF_EEGResponseFunctionGroup(sa, saveName)

groupDat = [];

if ~exist('saveName', 'var')
    saveName = input('name the output file', 's');
end

for s = 2:length(sa)
    %EF_MakeBetaMaps(sa{s})
    [S{s} res{s}] = EF_EEG_ResponseFunction(sa{s});
end

% for i = 1:length(res)
%     for t = 1:length(res{i}.classifyBOLDWithEEG.TR)
% 
%     groupDat.BOLDBP(i,t,:) = res{i}.tTestBOLDButtonPressVsFix.TR(t).stats.tstat;
%     groupDat.EEGBP(i,:) = res{i}.tTestEEGButtonPressVsFix.stats.tstat;
%     
%     groupDat.rRTWithEEG(i,:) = res{i}.classifyRTWithEEG.rObsVsPred;
%     groupDat.rBOLDWithEEG(i,t,:) = res{i}.classifyBOLDWithEEG.TR(t).rObsVsPred;
%     
%     groupDat.AccButtonPressWithEEG(i,:) = res{i}.classifyButtonPressWithEEG.meanAccuracyOverClasses;
%     
%     groupDat.pRTWithEEG(i,:) = res{i}.classifyRTWithEEG.pObsVsPred;
%     groupDat.pBOLDWithEEG(i,t,:) = res{i}.classifyBOLDWithEEG.TR(t).pObsVsPred;
%     end
% end
% 
subs = sa;
classMatFile = fullfile(S{1}.classMatDir, saveName);
save (classMatFile, 'S', 'res', 'subs');
%save (classMatFile, 'groupDat', 'S', 'res
% 
% for i=1:27, w1(i)=y.res{1}.classifyMemoryStateWithEEG; end
% for i=1:27, w2(i)=y.res{1}.classifyHCMemoryStateWithEEG; end
% for i=1:27, w3(i)=y.res{1}.classifyHCHitsVsHCCRs; end
% for i=1:27, w4(i)=y.res{1}.classifyHCHitsVsRecollection; end');
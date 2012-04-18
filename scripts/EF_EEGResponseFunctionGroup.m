function [groupDat res S] = EF_EEGResponseFunctionGroup(sa, saveName)

if ~exist('saveName', 'var')
    saveName = input('name the output file', 's');
end

for s = 1:length(sa)
    [S{s} res{s}] = EF_EEG_ResponseFunction(sa{s});
end

for i = 1:length(res)
    for t = 1:length(res{i}.classifyBOLDWithEEG.TR)

    groupDat.BOLDBP(i,t,:) = res{i}.tTestBOLDButtonPressVsFix.TR(t).stats.tstat;
    groupDat.EEGBP(i,:) = res{i}.tTestEEGButtonPressVsFix.stats.tstat;
    
    groupDat.rRTWithEEG(i,:) = res{i}.classifyRTWithEEG.rObsVsPred;
    groupDat.rBOLDWithEEG(i,t,:) = res{i}.classifyBOLDWithEEG.TR(t).rObsVsPred;
    
    groupDat.AccButtonPressWithEEG(i,:) = res{i}.classifyButtonPressWithEEG.meanAccuracyOverClasses;
    
    groupDat.pRTWithEEG(i,:) = res{i}.classifyRTWithEEG.pObsVsPred;
    groupDat.pBOLDWithEEG(i,t,:) = res{i}.classifyBOLDWithEEG.TR(t).pObsVsPred;
    end
end
% 
classMatFile = fullfile(S{1}.classMatDir, saveName);
save (classMatFile, 'groupDat', 'S', 'res');
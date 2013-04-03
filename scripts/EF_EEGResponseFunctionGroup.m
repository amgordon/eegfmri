function [groupDat res S] = EF_EEGResponseFunctionGroup(saveName, sa, varargin)


groupDat = [];

if ~exist('sa', 'var')
    sa = {'ef_091211' 'ef_091511' 'ef_092111' 'ef_092211' ...
        'ef_092711' 'ef_092911' 'ef_100511' 'ef_101411'...
        'ef_040512' 'ef_040712' 'ef_041112' 'ef_042912'};
elseif isempty(sa)
    sa = {'ef_091211' 'ef_091511' 'ef_092111' 'ef_092211' ...
        'ef_092711' 'ef_092911' 'ef_100511' 'ef_101411'...
        'ef_040512' 'ef_040712' 'ef_041112' 'ef_042912'};
end

if ~exist('saveName', 'var')
    saveName = input('What is the name of the output file? ', 's');
end

subs = sa;


for s = 1:length(sa)
    % EF_MakeBetaMaps(sa{s})
    [S{s} res{s}] = EF_EEG_ResponseFunction(sa{s}, varargin);
    classMatFile = fullfile(S{1}.classMatDir, saveName);
    save (classMatFile, 'S', 'res', 'subs', '-v7.3');
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


%fprintf('break');
%save (classMatFile, 'groupDat', 'S', 'res')
% 
% for i=1:27, w1(i)=y.res{1}.classifyMemoryStateWithEEG; end
% for i=1:27, w2(i)=y.res{1}.classifyHCMemoryStateWithEEG; end
% for i=1:27, w3(i)=y.res{1}.classifyHCHitsVsHCCRs; end
% for i=1:27, w4(i)=y.res{1}.classifyHCHitsVsRecollection; end');
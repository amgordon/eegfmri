function [] = EF_concatenateClassifierResults(filePattern)
cd /biac4/wagner/biac3/wagner5/alan/eegfmri/fmri_data/EEG_to_BOLD/classMats
d = dir([filePattern '*']);
dF = {d.name};

if length(dF)>17
   error('too many files match the filePattern'); 
end

for i = 1:length(dF)
    thisRes = load(dF{i});
    res.subj{i} = [thisRes.res{1}];
    %res.subjArray(i) = thisRes.res.subjArray;
end

save (filePattern, 'res');

for i = 1:length(dF)
    delete(dF{i})
end

end


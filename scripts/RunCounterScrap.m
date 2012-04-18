exp_dir = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/';

d = dir(fullfile(exp_dir, 'ef*'));

for i =1:length(d)
   cd (fullfile(exp_dir, (d(i).name), 'functional'))
   dS = dir('scan*');
   
   for s = 1:length(dS)
       cd (dS(s).name)
       maxScan_h = dir('scan*.nii');
       
       maxScan.sub(i).scan(s) = str2num(maxScan_h(end).name(9:11));
       cd ../
   end
   maxScan.sub(i).subName = d(i).name;
end
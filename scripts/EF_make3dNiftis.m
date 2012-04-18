function [] = EF_makeAnas(par)

d = dir(fullfile(par.rawdir, 'scan*'));


for i = 1:length(d)
    
    % make a run-specific diretory
    scanDir = ['scan' prepend(num2str(i))];
     newDir = fullfile(par.funcdir, scanDir);
     mkdir(newDir);
%     
%     % move the 4d nii to the run-specific directory
     movefile(fullfile(par.rawdir, d(i).name), newDir);
     cd (newDir);
%     
%     % 4d nii to 3d nii
     unix(['fslsplit ' d(i).name ' ' d(i).name(1:end-4) '_' ' -t']);
%     
%     % move raw scan back to its directory
     movefile(fullfile(par.funcdir, scanDir, d(i).name), par.rawdir);
%     
%     % unzip the newly made niis
     unix('gunzip scan*.gz');
    
    % hide discarded volumes
    discardDir = fullfile(par.funcdir, scanDir, 'discardedVols');
    mkdir(discardDir);
    dS = dir('scan*.nii');
    dropVols = {dS(1:par.dropvol).name}';
    
    for j = 1:length(dropVols)
        movefile(dropVols{j}, discardDir);
    end
end
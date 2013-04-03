function []=JR_make_mvpa_rawfilenames_list(subject)

par.allscanfiles = '';

topdir = pwd; % run from /Users/Jesse/fMRI/data/PAST/fMRI

par = par_params(subject);

% par.scans_to_include = [1 2 3 4 5 6 7 8 9 10];
% par.maxvol = 208; %highest volume number
% par.dropvol = 5; %dropped volumes (Assumed at start of scan)
% par.minvol = par.dropvol+1; %first non-dropped volume
% par.numvols = par.maxvol-par.dropvol; %number of volumes per scan
% par.funcdir = [topdir '/' subject '/functional'];  % or ..._functional_new, as appropriate

ac = 1;

for X = 1:length(par.scans_to_include)
    I = par.scans_to_include(X);
    for J = par.minvol:par.maxvol
        par.all_ra_scanfiles{ac} = [par.funcdir '/' prepend(num2str(I),2) '_func/ra' prepend(num2str(I),2) '_func.V' prepend(num2str(J),3) '.img']; %scan files
        par.all_sra_scanfiles{ac} = [par.funcdir '/' prepend(num2str(I),2) '_func/sra' prepend(num2str(I),2) '_func.V' prepend(num2str(J),3) '.img']; %ascan files
        par.all_sraf_scanfiles{ac} = [par.funcdir '/' prepend(num2str(I),2) '_func/sraf' prepend(num2str(I),2) '_func.V' prepend(num2str(J),3) '.img']; %wascan files
        par.all_wa_scanfiles{ac} = [par.funcdir '/' prepend(num2str(I),2) '_func/wa' prepend(num2str(I),2) '_func.V' prepend(num2str(J),3) '.img']; %swascan files
        par.all_s4mm_wa_scanfiles{ac} = [par.funcdir '/' prepend(num2str(I),2) '_func/s4mm_wa' prepend(num2str(I),2) '_func.V' prepend(num2str(J),3) '.img']; %swascan files
        par.all_s8mm_wa_scanfiles{ac} = [par.funcdir '/' prepend(num2str(I),2) '_func/s8mm_wa' prepend(num2str(I),2) '_func.V' prepend(num2str(J),3) '.img']; %swascan files
        ac = ac + 1;
    end
end

if ~exist([topdir '/' subject '/mvpa'],'dir')
    mkdir([topdir '/' subject '/mvpa'])
end

cd([topdir '/' subject '/mvpa'])

raw_filenames = par.all_ra_scanfiles;
save raw_filenames_ra.mat raw_filenames;

raw_filenames = par.all_sra_scanfiles;
save raw_filenames_sra.mat raw_filenames;

raw_filenames = par.all_sraf_scanfiles;
save raw_filenames_sraf.mat raw_filenames;

raw_filenames = par.all_wa_scanfiles;
save raw_filenames_wa.mat raw_filenames;

raw_filenames = par.all_s4mm_wa_scanfiles;
save raw_filenames_s4mm_wa.mat raw_filenames;

raw_filenames = par.all_s8mm_wa_scanfiles;
save raw_filenames_s8mm_wa.mat raw_filenames;

cd(topdir);
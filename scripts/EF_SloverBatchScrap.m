exp_dir = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data';
print_dir = '/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/slover_figures';


for i = 1:length(sa_good)
    
    maps = dir( (fullfile(exp_dir, sa_good{i}, 'analysis_hitsAndCRsByConf_Rec', 'mask.img')));
    
    anatFile = fullfile(exp_dir, sa_good{i}, 'anat', 'wV001.nii');
    InFile = fullfile(exp_dir, sa_good{i}, 'anat', 'In001.nii');
    
    if ~exist(anatFile)
        fname = InFile;
    else
        fname = anatFile;
    end
    %so.img(1).vol = spm_vol(fname);
    so.img(1).vol = spm_vol('/Applications/spm5/canonical/single_subj_T1.nii');
    
    so.img(2).vol = spm_vol(fullfile(exp_dir, sa_good{i}, 'analysis_hitsAndCRsByConf_Rec', 'mask.img'));
    %so.img(2).vol = spm_vol(fullfile(exp_dir, sa_good{i}, 'analysis_buttonPresses_byAmp_RTLocked', maps(1).name));
    %so.img(3).vol = spm_vol(fullfile(exp_dir, sa_good{i}, 'analysis_buttonPresses_byAmp_RTLocked', maps(1).name));
   

    figName = [sa_good{i}  'analysis_hitsAndCRsByConf_Rec_mask.eps'];
    slover(so);
    %so.img(3).range = -so.img(2).range;
    %so.img(1).range(2) = 4000;
    so.img(1).range(2) = 1;
    %%
    V = so.img(1).vol;
    D = V.dim(1:3);
    T = so.transform * V.mat;
    vcorners = [1 1 1; D(1) 1 1; 1 D(2) 1; D(1:2) 1; ...
        1 1 D(3); D(1) 1 D(3); 1 D(2:3) ; D(1:3)]';
    corners = T * [vcorners; ones(1,8)];
    SC = sort(corners');
    vxsz = sqrt(sum(T(1:3,1:3).^2));
    
    
    so.slicedef = [SC(1,1) vxsz(1) SC(8,1);SC(1,2) vxsz(2) SC(8,2)];
    
    so.slices = [SC(1,3):vxsz(3):SC(8,3)];
    stepSize = (SC(8,3) - SC(1,3) - 60)/20;
    so.slices = [(SC(1,3)+30):stepSize:(SC(8,3)-30)];
    %%
    
    so = paint(so);

    h = gcf;
    set(h, 'PaperPositionMode', 'auto');
    set(h, 'Position', [100 100 750 750]);
    print(h, '-depsc', fullfile(print_dir, figName));
    
end



% 
% base_dir = '/Users/alangordon/mounts/w5/alan/perceptMnemonic/fmri_data/group_analyses/patLinRegAnalysis_mnemonic_FaceHouseSeparate_SignedConf_14Subs';
% 
% f{1} = fullfile(base_dir, 'percFace_conf/spmT_0001.img');
% f{2} = fullfile(base_dir, 'percHouse_conf/spmT_0001.img');
% 
% dirChar = char(f{:});
% 
% 
% outputDir = fullfile(base_dir, 'greaterUnsignedEv_p01Conj');
% outputFile = fullfile(outputDir, 'spmT_0001.img');
% 
% if ~exist(outputDir)
%    mkdir(outputDir); 
% end
% 
% spm_imcalc_ui(dirChar,outputFile,'(i1>2.65) .* (i2<-2.65)');


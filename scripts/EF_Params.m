function par = EF_Params(substr)
% function par = parParams(subject)
%sets up parameters for parmap batching.  Currently set-up to run once per
%subject...

%profile on

par.FMCorrect = 1;
thisMachine = 'alan';
par.subTask = 'DM';
par.eegAnalysis = 1;
%% subject-specific stuff
switch substr
     case 'ef_071411'
        par.scansSelect = 1:3;
        par.goodEEGVols = par.scansSelect;
        par.numvols = [179 227 227];
        par.goodSub = 0;
        par.flagIt = 0;
        par.subNo = 1;
        par.sliceorder = [1:2:35 2:2:34]; % ventral to dorsal
    case 'ef_072111'
        par.scansSelect = [1 2 4 5];
        par.goodEEGVols = [1 2];
        par.numvols = [178 178 0 178 178];
        par.goodSub = 1;
        par.flagIt = 0;
        par.subNo = 2;
        par.sliceorder = [1:2:35 2:2:34]; % ventral to dorsal.
    case 'ef_083111'
        par.scansSelect = [1 4 5];
        par.goodEEGVols = par.scansSelect;
        par.numvols = [184 0 0 184 184];
        par.goodSub = 0;
        par.flagIt = 0;
        par.subNo = 3;
        par.sliceorder = [35:-2:1 34:-2:2];
    case 'ef_091211'
        par.scansSelect = 1:5;
        par.goodEEGVols = par.scansSelect;
        par.numvols = [234 234 234 234 234];
        par.goodSub = 1;
        par.flagIt = 0;
        par.subNo = 4;
        par.sliceorder = [35:-2:1 34:-2:2];
        par.criticalERPSamples = [500: 575];
    case 'ef_091311'
        par.scansSelect = 1:5;
        par.goodEEGVols = par.scansSelect;
        par.numvols = [225 225 225 225 225];
        par.goodSub = 0;
        par.flagIt = 0;
        par.subNo = 5;
        par.sliceorder = [35:-2:1 34:-2:2];
        par.criticalERPSamples = [];
    case 'ef_091511'
        par.scansSelect = 1:5;
        par.goodEEGVols = par.scansSelect;
        par.numvols = [225 225 225 225 225];
        par.goodSub = 1;
        par.flagIt = 0;
        par.subNo = 6;
        par.sliceorder = [35:-2:1 34:-2:2];
        par.criticalERPSamples = [500: 575];
    case 'ef_091911'
        par.scansSelect = 1:5;
        par.goodEEGVols = par.scansSelect;
        par.numvols = [225 225 225 225 234];
        par.goodSub = 0;
        par.flagIt = 0;
        par.subNo = 7;
        par.sliceorder = [35:-2:1 34:-2:2];
        par.criticalERPSamples = [];
    case 'ef_092111'
        par.scansSelect = 1:5;
        par.goodEEGVols = par.scansSelect;
        par.numvols = [225 225 225 225 225];
        par.goodSub = 1;
        par.flagIt = 0;
        par.subNo = 8;
        par.sliceorder = [35:-2:1 34:-2:2];
        par.criticalERPSamples = [500: 575];
    case  'ef_092211'
        par.scansSelect = 1:5;
        par.goodEEGVols = par.scansSelect;
        par.numvols = [225 225 225 225 225];
        par.goodSub = 1;
        par.flagIt = 0;
        par.subNo = 9;
        par.sliceorder = [35:-2:1 34:-2:2];
        par.criticalERPSamples = [463: 538];
    case  'ef_092711'
        par.scansSelect = [1 2 3 4 5];
        par.goodEEGVols = [1 2];
        par.numvols = [225 225 225 225 225];
        par.goodSub = 1;
        par.flagIt = 0;
        par.subNo = 10;
        par.sliceorder = [35:-2:1 34:-2:2];
        par.criticalERPSamples = [550: 625];
    case  'ef_092911'
        par.scansSelect = 1:5;
        par.goodEEGVols = par.scansSelect;
        par.numvols = [225 225 225 225 225];
        par.goodSub = 1;
        par.flagIt = 0;
        par.subNo = 11;
        par.sliceorder = [35:-2:1 34:-2:2];
        par.criticalERPSamples = [488: 563];
    case  'ef_100511'
        par.scansSelect = 1:5;
        par.goodEEGVols = par.scansSelect;
        par.numvols = [225 225 225 225 225];
        par.goodSub = 1;
        par.flagIt = 0;
        par.subNo = 6;
        par.subNo = 12;
        par.sliceorder = [35:-2:1 34:-2:2];
        par.criticalERPSamples = [463: 538];
    case 'ef_101411'
        %normalization is not good for this subject...
        % an effect of artefacts due to the eeg cap?
        par.scansSelect = 1:5;
        par.goodEEGVols = par.scansSelect;
        par.numvols = [225 225 225 225 225];
        par.goodSub = 1;
        par.flagIt = 0;
        par.subNo = 6;
        par.subNo = 13;
        par.sliceorder = [35:-2:1 34:-2:2];
        par.criticalERPSamples = [450: 525];
end

par.scans_to_include = par.scansSelect;

if (par.eegAnalysis)
    par.usedVols = par.numvols(par.goodEEGVols);
    par.usedScans = par.goodEEGVols;
else
    par.usedVols = par.numvols(par.scans_to_include);
    par.usedScans = par.scansSelect;
end

par.substr = substr;
%%
%----specify params-----:

% ---------Directory names----------
if strcmp(thisMachine, 'alan')
    par.exptdir = '/Users/alangordon/mounts/w5/alan/eegfmri/';
elseif strcmp(thisMachine, 'jesse')
    par.exptdir = '/Users/Jesse/w5/biac3/wagner5/alan/eegfmri/';
end

par.scriptsdir = fullfile(par.exptdir,'scripts');
par.fmridir = fullfile(par.exptdir, 'fmri_data');

par.subdir = fullfile(par.fmridir, par.substr);
par.anatdir = fullfile(par.subdir, 'anat');
par.funcdir = fullfile(par.subdir, 'functional');
par.behavdir = fullfile(par.subdir, 'behav');
par.rawdir = fullfile(par.subdir, 'raw');
par.logdir = fullfile(par.subdir, 'logfiles');
par.artrepdir = fullfile(par.subdir, 'artRepair');
par.classmat = fullfile(par.exptdir, '/fmri_data/mvpa_files/trainLocTestPerc_3Subs.mat');
par.meanfuncdir = fullfile(par.funcdir, 'meanFuncs');
par.analysisdir = fullfile(par.subdir, 'analysis_hitsVsCRs_HC');
par.ONAnalysisType = 'hitsVsCRs_HC';
par.erpFile = fullfile(par.subdir, 'erpData', 'trialdata_evonsetlock_LP125Hz.mat');
par.rawERPFile = fullfile(par.subdir, 'erpData', 'trialdata_RTlock.mat');
par.timeIntervalToInclude = 350:450;

%% exceptions


%%
% anat info
%dHiRes = dir(fullfile(par.anatdir, '*FSPGR*'));
%dInplane = dir(fullfile(par.anatdir, '*Inplane*'));

%par.hiresDir = fullfile(par.anatdir, dHiRes.name); 
%par.inplaneDir = fullfile(par.anatdir, dInplane.name); 
par.inplaneimg = [par.anatdir '/' 'In001.nii'];
par.hiresimg = [par.anatdir '/' 'V001.nii'];
%can add more later


%----------Scan Params----------------
%par.scans_to_include = [1 2 11];


par.numscans = sum(par.usedVols); %number of scans
par.art.motThresh = .5; % motion artefact threshold
par.art.sigThresh = 4; % signal intensity artefact threshold
par.TR = 2; %TR in seconds

par.numslice = 35;
par.TA = par.TR-par.TR/par.numslice; %TA
par.refslice = floor(par.numslice/2); %reference slice (assumes middle)

%par.maxvol = [126 600 600 600 600]; %highest volume number FOR EACH RUN
par.dropvol = 6; %dropped volumes (Assumed at start of scan)
par.minvol = par.dropvol+1; %first non-dropped volume
%par.numvols = par.maxvol-par.dropvol; %number of volumes per scan
par.slicetiming(1) = par.TA / (par.numslice -1);%timing var for slice timing
par.slicetiming(2) = par.TR - par.TA; %timing var for slice timing


% variables for realigning/reslicing

% realign flags...
par.realflag.quality = 0.9;
par.realflag.fwhm = 5;
par.realflag.sep = 4;
par.realflag.rtm = 1;
% par.realflag.PW;  %if field exists, then weighting is done...
par.realflag.interp = 2;

% reslice flags
par.reslflag.mask = 1;
par.reslflag.mean = 1;
par.reslflag.interp = 4;
par.reslflag.which = 2;

% re-orienting info
par.reorientmat = [105.0000 105.0000 59.4000 3.1416 0 -1.5708 1.0000 1.0000 1.0000 0 0 0]; %the proper rotation and translation to get the functional image preprocessed by kendrick's code to match the anatomical
dRN = dir([par.funcdir '/R*.nii']);
rawNiis_h = {};
for i = 1:length(dRN), rawNiis_h{i}=fullfile(par.funcdir, dRN(i).name); end %where the preprocessed .niis are stored
par.rawNiis = char(rawNiis_h{:});

% coregistration info
par.cor_meanfunc = [par.funcdir '/scan01/meanascan01_0006.nii']; %can change this to make more robust!
par.cor_inimg = par.inplaneimg;
par.cor_hiresimg = par.hiresimg;


% segment info
par.img2bSeg = par.cor_hiresimg;
par.segopts = struct('biascor',1,'GM',[0 0 1],'WM',[0 0 1],'CSF',[0 0 0],'cleanup',0);
par.segbiasfwhm = 60; % 60 is default in gui, 75 is default in command line for reasons unexplained
% see spm_config_preproc and spm_preproc(_write) for details


% normalization:
% gray matter:
par.graytemp = '/Applications/spm5/apriori/grey.nii';
par.grsegs(1,:) = [par.anatdir '/c1' 'V001.nii']; 
par.grsegs(2,:) = [par.anatdir '/c2' 'V001.nii']; 
par.graysrcimg = [par.anatdir '/c1' 'V001.nii']; 
par.graywrimg = [par.anatdir '/c1' 'V001.nii']; 
par.grflags.smosrc = 8;
par.grflags.smoref = 0;
par.grflags.regtype = 'mni';
par.grflags.cutoff = 25;
par.grflags.nits = 16;
par.grflags.reg = 1;
par.grwrflags.preserve = 0;
par.grwrflags.bb = [[-78 -112 -50];[78 76 85]]; % where did this come from?
par.grwrflags.vox        = [3 3 3];
par.grwrflags.interp     = 1;
par.grwrflags.wrap       = [0 0 0];

% spgm normalization:
par.spgrnrmflags.preserve = 0;
par.spgrnrmflags.bb = [[-78 -112 -50];[78 76 85]];
par.spgrnrmflags.vox = [2 2 2];
par.spgrnrmflags.interp     = 1;
par.spgrnrmflags.wrap       = [0 0 0];

% smoothing funcs
par.s4mm_smoothkernel = [4 4 4];
par.s8mm_smoothkernel = [8 8 8];
par.s6mm_smoothkernel = [6 6 6];

% specmaskvars

% specmaskvars
par.specwrflags.preserve = 0;
par.specwrflags.bb = [[-78 -112 -50];[78 76 85]];
par.specwrflags.vox        = [1 1 1];
par.specwrflags.interp     = 1;
par.specwrflags.wrap       = [0 0 0];

par.specsmooth = [8 8 8];%changed from [20 20 20] on 071409;
par.maskimg = [par.anatdir '/mask.nii']; 
par.smaskimg = [par.anatdir '/smask.nii'];
par.tsmaskimg = [par.anatdir '/tsmask.nii'];
par.wmaskimg = [par.anatdir '/wmask.nii'];
par.swmaskimg = [par.anatdir '/swmask.nii'];
par.tswmaskimg = [par.anatdir '/tswmask.nii'];
par.twmaskimg = [par.anatdir '/twmask.nii'];
par.addsegs = 'i1 + i2';
par.maskthresh = 'i1 > .2';


% model vars...
par.timing.fmri_t0 = par.refslice;  %micro-time resolution stuff (changed from 8)
par.timing.fmri_t = par.numslice; %used to be 16-- changed based on conversation with melina
par.timing.units = 'secs';
par.bases.hrf.derivs = [0 0]; % Melina says no cost to doing time and dispersion derivs ([0 0] = no derivs) NOTE that this will impact how you make contrasts!!
par.volt = 1;%if you don't want to model volterra interactions (and you don't want to model volterra interactions) leave this at 1.
% par.sess.scans is specified after populating list names...
par.sess.multi = {fullfile(par.behavdir, 'ons.mat')};
par.sess.multi_reg = {fullfile(par.behavdir, 'regs.mat')};
par.sess.hpf = 128;  % has anyone played around with this AND linear regs?
par.cvi = 'AR(1)'; %note that this actually gets changed to AR(0.2) in spm_fmri_spm_ui.  
% It looks like you might be able to give it a custom value
% by simply putting a number in here, but I haven't played around with it
par.global = 'None';
% explicit mask
par.mask = {};%{par.tswmaskimg};

% contrast vars
par.constat = 'T';

%par.srascanfiles = findScans(par, 'sRrun*.mat');
%par.rascanfiles = findScans(par, 'Rrun*.nii');
par.scanfiles  = findScans(par, 'scan*.nii');
par.ascanfiles  = findScans(par, 'ascan*.nii');
par.rascanfiles  = findScans(par, 'rascan*.nii');
par.srscanfiles  = findScans(par, 'srscan*.nii');
par.srascanfiles  = findScans(par, 'srascan*.nii');
par.rascanfilesPlusValidAndMean = findScans(par, 'rscan*');
par.wascanfiles.all = findScans(par, 'wascan*');
par.wrascanfilesPlusValidAndMean.all = findScans(par, 'rascan*.nii');
par.swascanfiles = findScans(par, 'swascan*.nii');

%spectral stuff
par.bandOfInterest = 1;
par.criticalSpectralSamples = 475:525; %95:155;
par.rawSpectralFile = fullfile(par.subdir, 'erpData', 'spectraldata3_RT.mat');
end

function scanfiles = findScans(par, matchstring)

scanfiles = [];

for i=par.usedScans;
    sf_h{i} = matchfiles (fullfile(par.funcdir, ['scan' prepend(num2str(i)) '/' matchstring ]));
end

sf_h2 = horzcat(sf_h{:});
scanfiles = vertcat(sf_h2{:});

end


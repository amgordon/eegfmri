function [pcregressors,r,params,bad,counts,snr,rrun] = fitprffast(stimulus,datafun,hrf,maxpolydeg,maxpcnum,wantresample,numinchunk,figuredir,xyzsize,omitparams)

% function [pcregressors,r,params,bad,counts,snr,rrun] = fitprffast(stimulus,datafun,hrf,maxpolydeg,maxpcnum,wantresample,numinchunk,figuredir,xyzsize,omitparams)
%
% <stimulus> is the stimulus with dimensions time x stim.  each column
%   should be zeros except for ones indicating trial onsets.  can be a cell 
%   vector whose elements represent different runs.  different runs can 
%   have different numbers of time points.
% <datafun> is a function that accept no arguments and returns either 
%   time x voxels or a cell vector of elements that are each time x voxels.
%   (different elements may have different numbers of time points.)
%   the data should be double or single format and should not contain any NaNs.
%   can also be the data itself (instead of a function).
% <hrf> is a column vector with the HRF.  the first point should be coicindent
%   with trial onset.  
% <maxpolydeg> is a non-negative integer with the maximum polynomial degree
%   to use for polynomial nuisance functions.
% <maxpcnum> is
%   N where N is the largest number of PCs to try.  in this case, we will 
%     evaluate PC numbers 0 through N, and the number of different fits we will obtain
%     is PCNUM = 1+N.
%  -A where A is the number of PCs to use.  in this case, the number of different
%     fits we will obtain is PCNUM = 2 (the first is the initial fit with 0 PCs).
%   X where X is a cell vector of matrices, each of which has PC regressors in 
%     the columns.  the same number of PC regressors should exist for each run.  
%     in this case, we will just use the PC regressors given by X as-is,
%     and the number of different fits we will obtain is PCNUM = 1.
% <wantresample> specifies the cross-validation scheme, and can be:
%   0 means no cross-validation (just fit fully)
%  -M where M is a matrix with dimensions A x B where A >= 1 indicates 
%     different cross-validation cases and B >= 2 matches the number of 
%     elements in the cell vector case of <stimulus>.  elements of M should
%     be consecutive integers, possibly skipping 0.  positive integers indicate runs to use 
%     for training; negative integers indicate runs to use for testing.  every cross-
%     validation case must involve at least one training run and one testing run.
%     for example, -[1 1 1 1 -1 -1 -2 -2 0] specifies a cross-validation scheme in which 
%     one cross-validation case is performed.  in this cross-validation case, the first four 
%     runs are used for training and then the 5th through 8th runs are used for testing.
%     additionally, we require that the concatenation of testing runs from the
%     cross-validation cases should exactly correspond to the full dataset.
% <numinchunk> is the number of voxels that you want to process simultaneously in the 
%   GLM fitting.  set this as high as possible, but according to how much memory you have.
%   some reasonable values might be 50000 or 100000.
% <figuredir> (optional) is the directory to write figures and results to.
%   if doesn't exist, we make it for you.  if [] or not supplied, do not 
%   write figures or results.
% <xyzsize> (optional) is a matrix size (e.g. [64 64 22]) that corresponds to the data.
%   this input matters only for the purposes of writing figures.  if this input is supplied, 
%   then we will attempt to interpret the data as a 3D matrix and we will write out 
%   various data inspections.  default: [].
% <omitparams> (optional) is whether to omit the <params> when saving the results (in order
%   to save disk space). default: 1.
%
% this function has some overlap in terms of functionality with fitprf.m.
% the big differences are that (1) this function can process a whole dataset in parallel
% instead of voxel by voxel and is thus much faster, (2) this function explicitly 
% implements the PCA-based de-noising strategy, and (3) this function is not as flexible
% as fitprf.m (i.e. it makes certain restrictions on the form of your data and model).
% 
% at a high-level, what we do is:
%   if <maxpcnum> is the N or -A case,
%     1. convolve the stimulus matrices with the HRF.
%     2. project out polynomial regressors from the stimulus and the data.
%     3. perform an initial GLM fit with no PCA regressors.
%     4. calculate the mean volume (i.e. the mean of each voxel's time-series),
%        and then determine which voxels are (1) bright (i.e. voxels that have a mean value
%        greater than 0.5 of the 99th percentile across all mean values) and (2) have no
%        stimulus-related modulation (i.e. voxels for which the initial GLM fit
%        resulted in 0 or less cross-validation R^2).
%     5. extract the time-series data for these voxels (after the polynomial projection),
%        normalize each time-series to be unit length, and then perform PCA, obtaining
%        the top <maxpcnum> (or abs(<maxpcnum>)) PCs.  there is one set of <maxpcnum> 
%        (or abs(<maxpcnum>)) PCs for each run.
%     6. re-fit the GLM model, either systematically varying the number of PCs from 1 to 
%        <maxpcnum> or setting the number of PCs to abs(<maxpcnum>).
%   if <maxpcnum> is the X case,
%     1. convolve the stimulus matrices with the HRF.
%     2. project out polynomial regressors from the stimulus and the data.
%     3. perform a GLM fit with the PCA regressors as given by <maxpcnum>.
%
% notes:
% - the R^2 metric is calculated using calccod.m (the <wantgain> and <wantmeansub> inputs
%   are set to 0).  note that model predictions reflect only the stimulus-related
%   portions of the model (not the polynomials nor the PCs).  also, the R^2 is computed
%   after polynomials are projected out from both the model prediction and the original data.
% - we take advantage of parfor, so turn on matlabpool if you can.
%
% upon completion, we return:
%  <pcregressors> as a cell vector of matrices that are each time x PCs
%  <r> as PCNUM x voxels with the R^2 values.  in the case that <maxpcnum> is the N case,
%    the first row reflects the initial GLM fit; the remaining rows reflect the 
%    PC-based fits from 1 through <maxpcnum>.  in the case that <maxpcnum> is the -A case,
%    the first row reflects the initial GLM fit; the second row reflects the PC-based fit.
%  <params> as parameters x voxels x resamples x PCNUM with the estimated GLM weights
%  <bad> as 1 x voxels with logical values indicating which voxels were selected for the PC analysis
%  <counts> as [A B C] where A is the number of voxels that were bright, B is the number of voxels
%    that had no stimulus-related modulation, and C is the number of voxels that satisfied both 
%    criteria (corresponding to the contents of <bad>).
%  <snr> as PCNUM x voxels with SNR values.
%    SNR is defined as:
%      snr = squish(max(abs(mean(params,3)),[],1) ./ sqrt(mean((std(params,[],3) * sqrt(size(params,3)-1)).^2,1)),3)';
%    this makes sense only if the resampling scheme is equivalent to jackknife resampling.
%    what we do is compute the standard error on the parameters (according to the jackknife),
%    then compute the pooled standard error, then compute the ratio between the maximum 
%    absolute parameter and this pooled standard error.  SNR values range from 0 to Inf.
%  <rrun> as PCNUM x voxels x runs with the R^2 values separated by runs.
% however, when <maxpcnum> is the X case, <bad> and <counts> will be returned as [].
%
% also, if <figuredir> is supplied, then we save the outputs of this function
% to a file named results.mat.  furthermore, we write out some useful figures:
%   1. summary.png shows the median R^2 across voxels for different number of PCs.
%   2. scatter01.png, scatter02.png, ... compare the R^2 obtained with no PCs
%      against the R^2 obtained with each different number of PCs.
%   3. meanvol.png shows the mean functional volume.
%   4. pcvoxels.png shows which voxels were used to derive PC regressors.
%   5. r00.png, r01.png, ... show the R^2 values for different numbers of PCs using
%      the range 0% to 100%.  (note that the colormap is nonlinearly scaled to 
%      enhance visibility.)
%   6. rsign00.png, ... show the R^2 values using the range -100% to 100%.
%   7. rthresh00.png, ... show which voxels have R^2 > 0.
%   8. in a directory called "runs", we write out
%      r00_run01.png, r00_run02.png, ..., r01_run01.png, r01_run02.png, ..., etc.
%      these show the R^2 values for different numbers of PCs and different runs using
%      the range 0% to 100%.  (note that the colormap is nonlinearly scaled to 
%      enhance visibility.)
%   9. snr00.png, snr01.png, ... show the SNR values for different numbers of PCs
%      using the range 0 to 10.
% note that when <maxpcnum> is the X case, the 'summary', 'scatter', 'meanvol',
% and 'pcvoxels' figures will not be produced.
%
% history:
% 2012/10/12 - implement speed-up; reduce memory usage.
% 2012/02/24 - allow <maxpcnum> to be -A.
% 2011/08/01 - allow <wantresample> to be 0; allow <maxpcnum> to be a cell vector
% 2011/07/26 - add <snr> output
% 2011/07/26 - add <rrun> output
% 2011/06/26 - add more figure output
% 2011/06/23 - add new inputs <figuredir>, <xyzsize>, <omitparams>
% 2011/04/10 - allow datafun to be the data itself; implement special mode for maxpcnum
% 2011/03/21 - first version

% INTERNAL NOTES:
% - maxpcnum can also be {pcregressors datasubset}.  this is a special mode for fitprffast_amountofdata.m.
%   pcregressors is cell vector of matrices, each of which has PC regressors in the columns.
%   datasubset is cell vector of indices (logical or 1-indices) indicating which data points
%   to use.  when we return, the only output that is valid is <r>, and it's a row vector.

% TODO: reduce memory usage, e.g. single
% TODO: we should x-val to select degree of polynomials!!

% internal constants
pthresh = 99;    % percentile for brain voxel values
pfrac = 0.5;     % what fraction of pthresh to set the threshold at

% input
if ~exist('figuredir','var') || isempty(figuredir)
  figuredir = [];
end
if ~exist('xyzsize','var') || isempty(xyzsize)
  xyzsize = [];
end
if ~exist('omitparams','var') || isempty(omitparams)
  omitparams = 1;
end
if ~iscell(stimulus)
  stimulus = {stimulus};
end
wantfigs = ~isempty(figuredir);

% calc
normalmode = ~iscell(maxpcnum);
weirdmode = iscell(maxpcnum) && length(maxpcnum)==2 && iscell(maxpcnum{2});
pcmode = iscell(maxpcnum) && ~weirdmode;

% make figure dir
if wantfigs
  mkdirquiet(figuredir);
  mkdirquiet([figuredir '/runs']);
end

% convolve stimulus matrices with HRF
fprintf('convolving stimulus matrices with HRF...'); stime = clock;
for p=1:length(stimulus)
  otime = size(stimulus{p},1);
  stimulus{p} = conv2(full(stimulus{p}),hrf);  % hrm, stimulus might be sparse, so convert here
  stimulus{p} = stimulus{p}(1:otime,:);
end
fprintf('done (%.1f minutes).\n',etime(clock,stime)/60);

% load in all data and convert to single [BIG MEMORY USAGE HERE]
fprintf('loading data...'); stime = clock;
if isa(datafun,'function_handle')
  data = feval(datafun);
else
  data = datafun;
end
if ~iscell(data)
  data = {data};
end
fprintf('done (%.1f minutes).\n',etime(clock,stime)/60);

% before we perturb the data, calculate mean volume for later
meanvol = meancell(data,1);  % 1 x voxels

% project out polynomials from both the data and stimulus.
% [can't use parfor; it complains about data being too large. but who cares because this loop is inherently multithreaded anyway.]
fprintf('projecting out polynomials from data and stimulus...'); stime = clock;
for p=1:length(data)
  pmatrix = constructpolynomialmatrix(size(data{p},1),0:maxpolydeg);
  if weirdmode
    pm = projectionmatrix(pmatrix(maxpcnum{2}{p},:));
    data{p} = pm*data{p}(maxpcnum{2}{p},:);
    stimulus{p} = pm*stimulus{p}(maxpcnum{2}{p},:);
  else
    pm = projectionmatrix(pmatrix);
    data{p} = pm*data{p};
    stimulus{p} = pm*stimulus{p};
  end
end
fprintf('done (%.1f minutes).\n',etime(clock,stime)/60);

% prep
dataval = [];  % this functionality is not being used yet....
dataclass = class(data{1});

% do it
if weirdmode

  % do the fitting
  fprintf('performing the GLM fit...'); stime = clock;
  r = fitglm(stimulus,data,cellfun(@(x,y) x(y,:),maxpcnum{1},maxpcnum{2},'UniformOutput',0),-size(maxpcnum{1}{1},2),wantresample,dataval,numinchunk);
  fprintf('done (%.1f minutes).\n',etime(clock,stime)/60);
  
  % deal with the other outputs
  pcregressors = [];
  params = [];
  bad = [];
  counts = [];
  snr = [];
  rrun = [];

elseif pcmode

  % do the fitting
  fprintf('performing the GLM fit...'); stime = clock;
  [r,params,rrun] = fitglm(stimulus,data,maxpcnum,-size(maxpcnum{1},2),wantresample,dataval,numinchunk);
  fprintf('done (%.1f minutes).\n',etime(clock,stime)/60);
  
  % deal with other outputs
  pcregressors = maxpcnum;
  bad = [];
  counts = [];

  % calculate snr
  snr = squish(max(abs(mean(params,3)),[],1) ./ sqrt(mean((std(params,[],3) * sqrt(size(params,3)-1)).^2,1)),3)';  % (1+maxpc) x voxels

else

  % do the initial fitting
  fprintf('performing initial GLM fit...'); stime = clock;
  [r,params,rrun] = fitglm(stimulus,data,[],0,wantresample,dataval,numinchunk);
  fprintf('done (%.1f minutes).\n',etime(clock,stime)/60);
  
  % figure out the bad voxels and calculate counts
  bad1 = meanvol > prctile(meanvol,pthresh)*pfrac;  % logical indicating voxels that are bright
  bad2 = r <= 0;                                    % logical indicating voxels that have no cross-validation power
  bad = bad1 & bad2;                                % logical indicating voxels that satisfy both criteria
  counts = [count(bad1) count(bad2) count(bad)];
  
  % do svd to determine pcregressors, which is a cell vector of matrices with PC regressors in the columns
  fprintf('calculating PC regressors...'); stime = clock;
  datatemp = cellfun(@(x) x(:,bad),data,'UniformOutput',0);  % subset this out because the parallelism was copying too much.  lame.
  pcregressors = {};
  parfor p=1:length(datatemp)
    pcregressors{p} = [];  % pre-declare to avoid stupid parfor complaints??
    [temp,len] = unitlengthfast(datatemp{p},1);
    temp = double(temp(:,len~=0));
    [u,s,v] = svds(temp*temp',abs(maxpcnum));
    pcregressors{p} = cast(u,dataclass);
  end
  clear datatemp temp len u s v;  % clear just to be sure
  fprintf('done (%.1f minutes).\n',etime(clock,stime)/60);
  
  % do the various fits
  if maxpcnum > 0
    fprintf('evaluating different numbers of PCs...'); stime = clock;
    [r(1+(1:maxpcnum),:),params(:,:,:,1+(1:maxpcnum)),rrun(1+(1:maxpcnum),:,:)] = ...
      fitglm(stimulus,data,pcregressors,maxpcnum,wantresample,dataval,numinchunk);
    fprintf('done (%.1f minutes).\n',etime(clock,stime)/60);
  elseif maxpcnum < 0
    fprintf('using PCs...'); stime = clock;
    [r(2,:),params(:,:,:,2),rrun(2,:,:)] = ...
      fitglm(stimulus,data,pcregressors,maxpcnum,wantresample,dataval,numinchunk);
    fprintf('done (%.1f minutes).\n',etime(clock,stime)/60);
  end

  % calculate snr
  snr = [];  % (1+maxpc) x voxels
  for pp=1:size(params,4)
    snr(pp,:) = flatten(max(abs(mean(params(:,:,:,pp),3)),[],1) ./ sqrt(mean((std(params(:,:,:,pp),[],3) * sqrt(size(params,3)-1)).^2,1)));
  end
  
end

% save results
if wantfigs
  if omitparams
    save([figuredir '/results.mat'],'pcregressors','r',         'bad','counts','snr','rrun');
  else
    save([figuredir '/results.mat'],'pcregressors','r','params','bad','counts','snr','rrun');
  end
end

% save figures
if wantfigs

  if ~pcmode
    % make figure showing summary
    pclist = choose(maxpcnum >= 0,1:maxpcnum-1,-maxpcnum);
    figureprep; hold on;
    bar(nanmedian(r,2)');
    set(gca,'XTick',1:size(r,1),'XTickLabel',[{'0'} mat2cellstr(pclist)]);
    xlabel('Number of PCs'); ylabel('Median R^2');
    title(sprintf('Counts = %s',mat2str(counts)));
    figurewrite('summary',[],[],figuredir);
  end
  
  if ~pcmode
    % make figure showing scatter plots
    pclist = choose(maxpcnum >= 0,1:maxpcnum,-maxpcnum);
    for p=1:length(pclist)
      figureprep([100 100 500 500]); hold on;
      scattersparse(r(1,:),r(1+p,:),50000,0,16,'r','.');
      axis([-50 100 -50 100]); axissquarify; axis([-50 100 -50 100]); 
      xlabel('Original R^2'); ylabel('New R^2');
      title(sprintf('Number of PCs = %d (showing at most 50,000 voxels)',pclist(p)));
      figurewrite(sprintf('scatter%02d',p),[],[],figuredir);
    end
  end
  
  % inspect pc voxels and r volumes
  if ~isempty(xyzsize)
    if ~pcmode
      imwrite(uint8(255*makeimagestack(reshape(meanvol,xyzsize),1)),gray(256),[figuredir '/meanvol.png']);
      imwrite(uint8(255*makeimagestack(reshape(bad,xyzsize),[0 1])),gray(256),[figuredir '/pcvoxels.png']);
    end
    for p=1:size(r,1)
      imwrite(uint8(255*makeimagestack(reshape(signedarraypower(r(p,:)/100,0.5),xyzsize),[0 1])),hot(256),sprintf([figuredir '/r%02d.png'],p-1));
      imwrite(uint8(255*makeimagestack(reshape(r(p,:),xyzsize),[-100 100])),cmapsign(256),sprintf([figuredir '/rsign%02d.png'],p-1));
      imwrite(uint8(255*makeimagestack(reshape(r(p,:)>0,xyzsize),[0 1])),gray(256),sprintf([figuredir '/rthresh%02d.png'],p-1));
      for q=1:size(rrun,3)
        imwrite(uint8(255*makeimagestack(reshape(signedarraypower(rrun(p,:,q)/100,0.5),xyzsize),[0 1])),hot(256),sprintf([figuredir '/runs/r%02d_run%02d.png'],p-1,q));
      end
    end
  end

  % inspect snr
  if ~isempty(xyzsize)
    for p=1:size(snr,1)
      imwrite(uint8(255*makeimagestack(reshape(snr(p,:),xyzsize),[0 10])),hot(256),sprintf([figuredir '/snr%02d.png'],p-1));
    end
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r,params,rrun] = fitglm(stimulus,data,pcregressors,pcnums,wantresample,dataval,numinchunk)

% function [r,params,rrun] = fitglm(stimulus,data,pcregressors,pcnums,wantresample,dataval,numinchunk)
%
% <stimulus> is a cell vector with the stimulus convolved with the HRF and with polynomials projected out
% <data> is a cell vector with all the data, with polynomials projected out.
% <pcregressors> is a cell vector of matrices.  each matrix has the PC regressors in the columns.
%   the number of columns must be equal to or greater than the max number in <pcnums>.
% <pcnums> is 0 or N or -N.  if 0, we don't use PCs (one fit, N=1).  
%   if N, we try PC numbers 1 through N (N fits).  if -N, we try PC number N (1 fit).
% <wantresample> is the original input to fitprffast.m
% <dataval> must be [] (NOT REALLY USED YET)
%
% for each number of PCs, fit the GLM model and quantify cross-validation R^2 (or full fit R^2)
% the R^2 values reflect stimulus-related activity but with polynomials projected out.
%
% we return:
%  <r> as N x voxels with the R^2 between the model predictions (or model fit)
%    and the data.  (polynomials are projected out from both
%    quantities before the R^2 calculation.)
%  <params> as parameters x voxels x resamples x N with
%    the weights on the different stimuli.
%  <rrun> as N x voxels x runs with the R^2 between the model 
%    predictions (or model fit) and the data, separated by runs.
%    (polynomials are projected out from both quantities before the 
%    R^2 calculation.)

% calc
nparams = size(stimulus{1},2);      % number of free parameters
nvoxels = size(data{1},2);          % number of voxels
npcs = choose(pcnums<=0,1,pcnums);  % how many different numbers of PCs to try
nresamples = size(wantresample,1);  % how many different resamplings are we going to do?
dataclass = class(data{1});         % what format, like 'single'

% initialize
r = zeros(npcs,nvoxels,dataclass);
rrun = zeros(npcs,nvoxels,length(data),dataclass);
params = zeros(nparams,nvoxels,nresamples,npcs,dataclass);  % WARNING: THIS TAKES UP QUITE A BIT OF MEMORY
pctodo = choose(pcnums==0,[0],choose(pcnums>0,1:pcnums,[-pcnums]));
chunks = chunking(1:nvoxels,numinchunk);

% do it [TODO: FOR A SPEEDUP, CAN'T WE PRECOMPUTE AND CACHE THE STIMULUS-RELATED MATRICES HERE?]
for cc=1:length(chunks)
  fprintf('*** processing chunk %d of %d ***\n',cc,length(chunks));

  % initialize second copy of data and stimulus [painful but necessary]
  dataB = cellfun(@(x) x(:,chunks{cc}),data,'UniformOutput',0);  % aha, dataB is already subseted
  stimulusB = stimulus;
  
  % do it
  for z=1:length(pctodo)
    fprintf('evaluating pc number=%d. ',pctodo(z));
  
    % prepare stimulus and data (which have the PCs projected out, if applicable)
    if pctodo(z) ~= 0
      fprintf('projecting out PCs from data and stimulus...'); stime = clock;
      for zz=1:length(stimulus)  % USE PARFOR??
        if pcnums > 0
          iii = pctodo(z);
        else
          iii = 1:pctodo(z);  % need to use the whole range in one shot
        end
        pcrmatrix = pcregressors{zz}(:,iii);
%%        pcrmatrix = pcrmatrix(:,~all(pcrmatrix==0,1));
        ppp = projectionmatrix(pcrmatrix);  % prepare projection matrix for the current PC
        stimulusB{zz} = ppp*stimulusB{zz};  % project it out from the stimulus
        dataB{zz} = ppp*dataB{zz};  % project it out from the data
      end
      fprintf('done (%.1f minutes).\n',etime(clock,stime)/60);
    end
  
    % fit model and calculate model predictions
    modelfit = {};  % cell vector of dimensions 1 x runs.  each element is time x voxels.  NOTE: heavy memory usage!
    if isequal(wantresample,0)
      fprintf('fit fully...'); stime = clock;
      trainix = logical(ones(1,length(data)));
      params(:,chunks{cc},1,z) = mtimescell(olsmatrix(cat(1,stimulusB{trainix})),dataB(trainix));
      modelfit = cellfun(@(x) x*params(:,chunks{cc},1,z),stimulus(trainix),'UniformOutput',0);  % note stimulus not stimulusB
      fprintf('done (%.1f minutes).\n',etime(clock,stime)/60);
    else
      for p=1:size(wantresample,1)  % DON'T PARFOR BECAUSE HEAVY MEMORY USAGE
        fprintf('fit and predict for resampling case %d...',p); stime = clock;
        trainix = wantresample(p,:) < 0;
        testix = wantresample(p,:) > 0;
        params(:,chunks{cc},p,z) = mtimescell(olsmatrix(cat(1,stimulusB{trainix})),dataB(trainix));
        modelfit = [modelfit cellfun(@(x) x*params(:,chunks{cc},p,z),stimulus(testix),'UniformOutput',0)];  % note stimulus not stimulusB
        fprintf('done (%.1f minutes).\n',etime(clock,stime)/60);
      end
    end

% % check shifts [EXPERIMENTAL]
% tryits=[-10:.5:10];  % or [-3:.5:3]
% for ppp=1:length(modelfit)  % for each run
%   for qqq=1:length(tryits)  % for each shift
%     temp = calccod(modelfit{ppp},sincshift(data{ppp},tryits(qqq),1),1,0,0);
%     imwrite(uint8(255*makeimagestack(reshape(signedarraypower(temp/100,0.5),[80 80 28]),[0 1])), ...
%       hot(256),sprintf('shift_run%02d_%02d.png',ppp,qqq));
%   end
% end
  
    % calculate R^2 goodness
    fprintf('calculating R^2...'); stime = clock;
    r(z,chunks{cc}) = calccodcell(modelfit,cellfun(@(x) x(:,chunks{cc}),data,'UniformOutput',0),1);  % note we use data not dataB here!
    rrun(z,chunks{cc},:) = catcell(3,cellfun(@(aa,bb) calccod(aa,bb,1,0,0),modelfit,cellfun(@(x) x(:,chunks{cc}),data,'UniformOutput',0),'UniformOutput',0));
              %%choose(isempty(dataval),data,dataval)
    
    clear modelfit;  % let's clean up
    fprintf('done (%.1f minutes).\n',etime(clock,stime)/60);
  
  end
  clear stimulusB dataB;  % finally, we can clean these up

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUNK, OBSOLETE STUFF

% % input
% if ~exist('smoothfwhm','var') || isempty(smoothfwhm)
%   smoothfwhm = [];
% end
% 
% % smooth data if requested    [THIS WILL BE INCOMPATIBLE WITH FLATTENED FORMAT, RE CHECK AND RE IMPLEMENT]
% if isempty(smoothfwhm)
% else
%   dataval = data;
%   data = smoothvolumes(data,[1 1 1],smoothfwhm);
% end
% 
%   if ~isempty(smoothfwhm)  % RECHECK
%     dataval{p} = pm*squish(dataval{p},3)';  % time x voxels
%   end
% 
% % <NONOmetric>.    r(z,:) = calccod(catcell(1,modelfit),catcell(1,data),1,0,0);
% %   else
% %     r(z,:) = calccod(catcell(1,modelfit),catcell(1,dataval),1,0,0);
% % <NONOmetric> (optional) is a function that accepts two column vectors (the first is the model,
% %   the second is the data) and outputs a value.  default is @calccod, which computes the 
% %   percent variance of the data explained by the model.  note that polynomial nuisance functions
% %   are projected out from both the data and the model before computation of the <metric>
% %   (for more details, see below).




% % make figuredir if necessary
% mkdirquiet(figuredir);
% figuredir = subscript(matchfiles(figuredir),1,1);




% % RESIDUAL MODELING
% % 
% %   % for each run
% [baddim1,baddim2,baddim3] = ind2sub(xyzsize,find(bad));
% slicedata = NaN*zeros(xyzsize(1:2));
% assert(xyzsize(1)==xyzsize(2));
% presid = {};
% for p=1:length(datatemp)
%   presid{p} = zeros(xyzsize(1),xyzsize(2),xyzsize(3),size(datatemp{p},1));
%   % for each time point
%   for q=1:size(datatemp{p},1)
%     % for each slice
%     for r=1:xyzsize(3)
%       slicedata(:) = NaN;
%       iii = sub2ind(xyzsize(1:2),baddim1(baddim3==r),baddim2(baddim3==r));
%       slicedata(iii) = datatemp{p}(q,baddim3==r);
%       
%       [xx,yy] = ndgrid(1:100,1:100);
%       whoa = [1.1 1.5 2 3 5 7 10];
%       ress = [];
%       for zzz=1:length(whoa), zzz
%         prediction = slicedata;
%         for ppp=1:50, ppp
%           [rr,d,rrnot] = picksubset(iii,[50 ppp]);
%           okok = localregression2d(xx(rrnot),yy(rrnot),slicedata(rrnot),xx(rr),yy(rr),[],[],whoa(zzz));
%           prediction(rr) = okok;
%         end
%         ress(zzz) = calccorrelation(prediction(iii),slicedata(iii));
%       end
%       figure;bar(ress);
% 
%       whoa = [1.1 1.3 1.5 1.7 2 2.5 3];
%       ress = [];
%       for zzz=1:length(whoa)
%         okok = localregression2d(xx(iii),yy(iii),slicedata(iii),xx,yy,[],[],whoa(zzz));
%         ress(zzz) = calccorrelation(okok(iii),slicedata(iii));
%       end
%       figure;bar(ress);
% 
%     end
%   end
% end
% 
% 
% 
%       figure; imagesc(slicedata,[-100 100]);
%       
% 
% grats = makegratings2d(xyzsize(1),[1 2 3 4 5 6],40,2);
% grats(:,:,end+1) = 1;  % res x res x regressors
% grats = squish(grats,2);
% 
%       hh = pinv(grats(iii,:))*vflatten(slicedata(iii));
%       modelfit = grats*hh;
%       
%       figure;
%       scatter(slicedata(iii),modelfit(iii),'r.');
%       
%       p

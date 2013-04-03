function results = fitprffast_amountofdata(stimulus,datafun,hrf,maxpolydeg,wantresample,numinchunk,pcregressors,datapct,numrep)

% function results = fitprffast_amountofdata(stimulus,datafun,hrf,maxpolydeg,wantresample,numinchunk,pcregressors,datapct,numrep)
%
% <stimulus>,<datafun>,<hrf>,<maxpolydeg>,<wantresample>,<numinchunk> are as in fitprffast.m
% <pcregressors> is a cell vector of matrices that are each time x PCs
% <datapct> (optional) is a vector of fractions to try.  each fraction should be in (0,1].
%   default: [0.2 0.4 0.6 0.8 1.0];
% <numrep> (optional) is the number of resamplings to perform for each fraction.  default: 2.
%
% this function is essentially a wrapper for fitprffast.m.  the purpose is to evaluate the effect
% of the amount of data on the cross-validation R^2 values.  
%
% for each fraction listed in <datapct>, we fit the GLM model <numrep> times.  each time, 
% we pick a random subset of the data points in each run (seeding based on clock.m) according
% to the fraction desired (rounding to the nearest whole number).  this involves
% subsetting the data, the stimulus, and the PC regressors.  then, in the fitting process, the
% cross-validation scheme specified by <wantresample> is carried out as usual.
%
% we return <results> as a matrix with dimensions length(<datapct>) x voxels.
% each row indicates the median R^2 values obtained for a given fraction of the data.
% (the median is taken across the <numrep> repetitions.)
%
% note that we do not evaluate different numbers of PCs (like in fitprffast.m);
% instead we always use the PCs specified in <pcregressors>.  also, note that
% we perform the subsetting before the cross-validation; thus, not only is there
% less data to fit the GLM model, but there is also less data to evaluate 
% cross-validation performance.
%
% history:
% 2011/04/10 - first version

% PERHAPS UPDATE THIS FOR NEW SNR OUTPUTS?

% input
if ~exist('datapct','var') || isempty(datapct)
  datapct = .2:.2:1;
end
if ~exist('numrep','var') || isempty(numrep)
  numrep = 2;
end
if ~iscell(stimulus)
  stimulus = {stimulus};
end

% load data once and for all [NOTE: MIRRORED FROM FITPRFFAST.M]
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

% do it
results = zeros(length(datapct),numrep,size(data{1},2));
for p=1:length(datapct)
  fprintf('***** processing datapct %d of %d (%.5f) *****\n',p,length(datapct),datapct(p));
  for q=1:numrep
    fprintf('******* rep %d of %d *****\n',q,numrep);
    
    % calculate datasubset
    datasubset = {};
    for r=1:length(stimulus)
      n = size(stimulus{r},1);
      datasubset{r} = picksubset(1:n,round(datapct(p)*n),sum(100*clock));
    end
    
    % do it
    [d,results(p,q,:)] = fitprffast(stimulus,data,hrf,maxpolydeg,{pcregressors datasubset},wantresample,numinchunk);

  end
end
results = squish(median(results,2),2);

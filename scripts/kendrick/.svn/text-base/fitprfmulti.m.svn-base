function fitprfmulti(opt,chunksize,chunknum)

% function fitprfmulti(opt,chunksize,chunknum)
%
% <opt> is a struct with the following fields:
%   <outputdir> is the directory to save results to
%   <fitprfinputs> is a cell vector with inputs to fitprf.m.  the <response> input
%     to fitprf.m should be passed here as [] (since the data are handled separately
%     by the <datafun> input).
%   <datafun> is a function that accepts a vector of voxel indices and returns
%     either voxels x time or a cell vector of elements that are each voxels x time.
%     (different elements may have different numbers of time points.)  the data
%     should be double-format.
%   <vxsfile> is a .mat file from which we get voxel indices 'vxs' from.
%     alternatively, can be a vector of voxel indices (row or column vector).
%     we automatically sort the indices and ensure uniqueness.
%     also, can be a .nii file --- in this case, we simply perform find() on
%     the volume to figure out what the voxel indices are.
%   <hrffile> (optional) is a .mat file from which we get 'hrfs' from.  this should be
%     a matrix t x length(vxs) with the hrf to assume for each voxel.  there should
%     be a 1-to-1 correspondence between the sorted, unique voxel indices specified
%     by <vxsfile> and the columns of 'hrfs'.  when <hrffile> is supplied,
%     the 'hrfmodel' parameter in <fitprfinputs> is ignored.  default is [] which
%     means to do the normal behavior (i.e. use the 'hrfmodel' parameter in <fitprfinputs>).
%   <wantsignaldrift> (optional) is whether to save the 'signal' and 'drift' outputs of fitprf.m.
%     we save these as single to save on space.  default: 0.
% <chunksize> is the number of voxels to do in a chunk
% <chunknum> is the chunk to process.  note that the maximum number of chunks is 999999 (since we use %06d).
%
% fit PRF models.  variables saved:
%  'params' is A x B x V
%  'paramsse' is 1 x A x V
%  'r' is 1 x V
%  'rrun' is A x V
%  'polyparams' is A x B x C x V
%  'polymeans' is A x B x V
%  'numiters' is [] or 1 x A x V
%  'hrf' is A x B x V
%  'betas' is A x B x V
%  'signal' is A x V
%  'drift' is A x V
% where V indicates different voxels and the other letters are placeholders for matrix dimensions
% (not necessarily the same values for different variables).
%
% history:
% 2011/06/23 - allow <vxsfile> to be a .nii file
% 2010/12/02 - the <meanint>,<driftstd>,<signalrms>,<noiserms> outputs have been removed.
%              remove <arfile> input (obsolete).
% 2010/10/11 - signal and drift saved as single.
% 2010/10/09 - add meanint,driftstd,signalrms,noiserms outputs.
% 2010/10/06 - add signal and drift outputs from fitprf.m and <wantsignaldrift>.
% 2010/09/30 - add reporting of time for each voxel's fit
% 2010/09/08 - update to be consistent with update of fitprf.m.
% 2010/06/15 - fitprf uses calccod now.  update to be consistent with that.
% 2010/05/06 - first version (modified from fitrf.m)

%%%%%%%%%%%%%%%% INPUTS

if ~isfield(opt,'hrffile')
  opt.hrffile = [];
end
if ~isfield(opt,'wantsignaldrift')
  opt.wantsignaldrift = 0;
end

%%%%%%%%%%%%%%%% PREPARATION

% make outputdir if necessary
if ~exist(opt.outputdir,'dir')
  mkdirquiet(opt.outputdir);
end
opt.outputdir = subscript(matchfiles(opt.outputdir),1,1);
outputfile = sprintf([opt.outputdir '/%06d.mat'],chunknum);

% report
fprintf(['*** fitprfmulti called:\n' ...
         '*** outputdir = %s\n' ...
         '*** chunksize = %d\n' ...
         '*** chunknum = %d\n'], ...
  opt.outputdir,chunksize,chunknum);
stime = clock;  % start time

% prepare voxels to process
if ischar(opt.vxsfile)
  if isequal(getextension(opt.vxsfile),'.mat')
    vxsfull = sort(union([],flatten(loadmulti(opt.vxsfile,'vxs'))));
  elseif isequal(getextension(opt.vxsfile),'.nii')
    temp = load_untouch_nii(opt.vxsfile);
    vxsfull = flatten(find(temp.img));
  else
    die;
  end
else
  vxsfull = sort(union([],flatten(opt.vxsfile)));
end
[vxs,vxbegin,vxend] = chunking(vxsfull,chunksize,chunknum);
vnum = length(vxs);

% deal with hrffile
if ~isempty(opt.hrffile)
  hrfs = loadmulti(opt.hrffile,'hrfs');
  hrfs = hrfs(:,vxbegin:vxend);
end

% save initial time
clear vxsfull;  % clear some big unnecessary variables
save(outputfile);

%%%%%%%%%%%%%%%% DO IT

% load in stimulus and extraregressors once and for all
if isa(opt.fitprfinputs{1},'function_handle')
  opt.fitprfinputs{1} = feval(opt.fitprfinputs{1});
end
if isa(opt.fitprfinputs{13},'function_handle')
  opt.fitprfinputs{13} = feval(opt.fitprfinputs{13});
end

% load data
data = feval(opt.datafun,vxs);

% loop
params = []; paramsse = []; r = []; rrun = []; polyparams = []; polymeans = []; 
numiters = []; hrf = []; betas = []; signal = []; drift = []; 
for p=1:vnum
  fprintf('now processing voxel %d (%d of %d).\n',vxs(p),p,vnum); vtime = clock;
  
  % prepare inputs to fitprf.m
  inputs = opt.fitprfinputs;
  if iscell(data)
    inputs{2} = cellfun(@(x) x(p,:)',data,'UniformOutput',0);
  else
    inputs{2} = data(p,:)';
  end
  if ~isempty(opt.hrffile)
    inputs{4} = hrfs(:,p);
  end

  % call fitprf.m
  if p==1  % we have to collect outputs like this to avoid MATLAB's automatic vector reshaping issue
    [params,paramsse,r,rrun,polyparams,polymeans,numiters,hrf,betas,signal,drift] = fitprf(inputs{:});
    rrun = rrun';  % make into column vector for outputing from this function
  else
    [params(:,:,p),paramsse(:,:,p),r(p),rrun(:,p),polyparams(:,:,:,p),polymeans(:,:,p), ...
     numiters0,hrf0,betas0,signal(:,p),drift(:,p)] = fitprf(inputs{:});

    % [] is a stupid special case that needs to be handled
    if ~isempty(numiters0)
      numiters(:,:,p) = numiters0;
    end
    if ~isempty(hrf0)
      hrf(:,:,p) = hrf0;
    end
    if ~isempty(betas0)
      betas(:,:,p) = betas0;
    end

  end
  
  fprintf('the fit for voxel %d (%d of %d) took %.1f seconds.\n',vxs(p),p,vnum,etime(clock,vtime));
end

% convert to save storage
signal = single(signal);
drift = single(drift);

% save
if opt.wantsignaldrift
  save(outputfile,'params','paramsse','r','rrun','polyparams','polymeans','numiters','hrf','betas','signal','drift','-append');
else
  save(outputfile,'params','paramsse','r','rrun','polyparams','polymeans','numiters','hrf','betas','-append');
end

%%%%%%%%%%%%%%%% DONE

fprintf('*** fitprfmulti completed in %.1f minutes.\n',etime(clock,stime)/60);









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUNK

% isbadcase = zeros(1,vnum);
%     isbadcase(p) = any(cellfun(@(x) any(isnan(x)),inputs{2}));
%     isbadcase(p) = any(isnan(inputs{2}));
%   if ~isbadcase(p)
% 
% % for bad cases, need to fill in NaNs in the results
% params(:,:,isbadcase) = NaN;
% paramsse(:,:,isbadcase) = NaN;
% r(isbadcase) = NaN;
% rrun(:,isbadcase) = NaN;
% polyparams(:,:,:,isbadcase) = NaN;
% polymeans(:,:,isbadcase) = NaN;
% numiters(:,:,isbadcase) = NaN;
% 
% ,'isbadcase
%
% REMOVED!
% %  'meanint' is A x V
% %  'driftstd' is A x V
% %  'signalrms' is A x V
% %  'noiserms' is A x V
% meanint = []; driftstd = []; signalrms = []; noiserms = [];
% ,meanint(:,p),driftstd(:,p),signalrms(:,p),noiserms(:,p)
% 
%  ...
%                   'meanint','driftstd','signalrms','noiserms',
%                                     ...
%                   'meanint','driftstd','signalrms','noiserms',




%   <arfile> (optional) is a .mat file from which we get 'ars' from.  this should be
%     a matrix 1 x length(vxs) with the AR value to assume for each voxel.  there should
%     be a 1-to-1 correspondence between the sorted, unique voxel indices specified
%     by <vxsfile> and the columns of 'ars'.  when <arfile> is supplied,
%     the 'ar' parameter in <fitprfinputs> is ignored.  default is [] which
%     means to do the normal behavior (i.e. use the 'ar' parameter in <fitprfinputs>).

% if ~isfield(opt,'arfile')
%   opt.arfile = [];
% end
% 
% % deal with arfile
% if ~isempty(opt.arfile)
%   ars = loadmulti(opt.arfile,'ars');
%   ars = ars(vxbegin:vxend);
% end
% 
%   if ~isempty(opt.arfile)
%     inputs{7} = ars(p);
%   end

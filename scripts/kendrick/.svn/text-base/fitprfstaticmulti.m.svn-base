function fitprfstaticmulti(opt,chunksize,chunknum)

% function fitprfstaticmulti(opt,chunksize,chunknum)
%
% <opt> is a struct with the following fields:
%   <outputdir> is the directory to save results to
%   <fitprfstaticinputs> is a cell vector with inputs to fitprfstatic.m.  the <response> input
%     to fitprfstatic.m should be passed here as [] (since the data are handled separately
%     by the <datafun> input).
%   <datafun> is a function that accepts a vector of voxel indices and returns
%     voxels x data x trials.
%   <vxsfile> is a .mat file from which we get voxel indices 'vxs' from.
%     alternatively, can be a vector of voxel indices (row or column vector).
%     we automatically sort the indices and ensure uniqueness.
%     also, can be a .nii file --- in this case, we simply perform find() on
%     the volume to figure out what the voxel indices are.
% <chunksize> is the number of voxels to do in a chunk
% <chunknum> is the chunk to process.  note that the maximum number of chunks is 999999 (since we use %06d).
%
% fit static PRF models.  variables saved:
%  'params' is A x B x V
%  'paramsse' is 1 x A x V
%  'r' is 1 x V
%  'modelfit' is A x V
%  'numiters' is [] or 1 x A x V
%  'responseB' is A x V
%  'rTRAIN' is 1 x V
%  'modelfitTRAIN' is C x V
%  'responseTRAIN' is C x V
% where V indicates different voxels and the other letters are placeholders for matrix dimensions
% (not necessarily the same values for different variables).
%
% history:
% 2011/01/23 - implement code speedup in a special case
% 2011/09/22 - in fitprfstatic.m, add <wanttrain> and the corresponding output variables
% 2011/08/01 - allow <vxsfile> to be a .nii file
% 2010/08/26 - implement some speed-ups
% 2010/08/12 - add reporting of time for each voxel's fit; add responseB
% 2010/08/10 - first version (modified from fitprfmulti.m)

%%%%%%%%%%%%%%%% INPUTS

%%%%%%%%%%%%%%%% PREPARATION

% make outputdir if necessary
if ~exist(opt.outputdir,'dir')
  mkdirquiet(opt.outputdir);
end
opt.outputdir = subscript(matchfiles(opt.outputdir),1,1);
outputfile = sprintf([opt.outputdir '/%06d.mat'],chunknum);

% report
fprintf(['*** fitprfstaticmulti called:\n' ...
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

% save initial time
save(outputfile);

%%%%%%%%%%%%%%%% DO IT

% load in stimulus once and for all
if isa(opt.fitprfstaticinputs{1},'function_handle')
  opt.fitprfstaticinputs{1} = feval(opt.fitprfstaticinputs{1});
end

% massage input so we can do the next step in one shot
if iscell(opt.fitprfstaticinputs{3}) && ~iscell(opt.fitprfstaticinputs{3}{1})
  opt.fitprfstaticinputs{3} = {opt.fitprfstaticinputs{3}};
end

% perform memory saving and code speedup in a very special case
if ~isequal(opt.fitprfstaticinputs{1},0) && ...  % if stimulus is not 0
   iscell(opt.fitprfstaticinputs{3}) && ...  % if model is the cell vector case
   length(opt.fitprfstaticinputs{3}{1}) >= 4 && ...  % if there is a W input
   ~isempty(opt.fitprfstaticinputs{3}{1}{4})  % if that W input is not empty
  test = feval(opt.fitprfstaticinputs{3}{1}{4},opt.fitprfstaticinputs{1});  % perform the transformation
  if ~iscell(test)  % if the transformation did not result in a cell
    fprintf('*** fitprfstaticmulti: performing stimulus transformation speed-up. ***\n');
    opt.fitprfstaticinputs{1} = test;   % ah, we can substitute in the transformed stimulus
    opt.fitprfstaticinputs{3}{1}{4} = [];  % we no longer need the transformation
  end
end

% load data
data = feval(opt.datafun,vxs);

% loop
params = []; paramsse = []; r = []; modelfit = []; numiters = []; responseB = [];
rTRAIN = []; modelfitTRAIN = []; responseTRAIN = [];
for p=1:vnum
  fprintf('now processing voxel %d (%d of %d).\n',vxs(p),p,vnum); vtime = clock;
  
  % prepare inputs to fitprf.m
  opt.fitprfstaticinputs{2} = squish(data(p,:,:),2);

  % call fitprfstatic.m
  if p==1  % we have to collect outputs like this to avoid MATLAB's automatic vector reshaping issue
    [params,paramsse,r,modelfit,numiters,responseB,rTRAIN,modelfitTRAIN,responseTRAIN] = fitprfstatic(opt.fitprfstaticinputs{:});
  else
    [params(:,:,p),paramsse(:,:,p),r(p),modelfit(:,p),numiters0,responseB(:,p),rTRAIN0,modelfitTRAIN0,responseTRAIN0] = fitprfstatic(opt.fitprfstaticinputs{:});
    if ~isempty(numiters0)  % [] is a special case that needs to be handled
      numiters(:,:,p) = numiters0;
    end
    if ~isempty(rTRAIN0)  % [] is a special case that needs to be handled
      rTRAIN(p) = rTRAIN0;
    end
    if ~isempty(modelfitTRAIN0)  % [] is a special case that needs to be handled
      modelfitTRAIN(:,p) = modelfitTRAIN0;
    end
    if ~isempty(responseTRAIN0)  % [] is a special case that needs to be handled
      responseTRAIN(:,p) = responseTRAIN0;
    end
  end
  
  fprintf('the fit for voxel %d (%d of %d) took %.1f seconds.\n',vxs(p),p,vnum,etime(clock,vtime));
end

% save
save(outputfile,'params','paramsse','r','modelfit','numiters','responseB','rTRAIN','modelfitTRAIN','responseTRAIN','-append');

%%%%%%%%%%%%%%%% DONE

fprintf('*** fitprfstaticmulti completed in %.1f minutes.\n',etime(clock,stime)/60);

function greconsmultifreq(grecons,files,figuredir,n,s,llbands,mode,numsteps,numframes,todo,dfactor,isinout,gain,fvol,dchack,tweaks,usesgerun)

% function greconsmultifreq(grecons,files,figuredir,n,s,llbands,mode,numsteps,numframes,todo,dfactor,isinout,gain,fvol,dchack,tweaks,usesgerun)
% 
% <grecons> is a string indicating the grecons utility that you wish to use.
%   the following versions have been tested: grecons* rev 95.
%   no guarantees if you use a different version of grecons.
% <files> is a file pattern (see matchfiles.m) matching one or more .7 files.
%   each .7 file should correspond to one functional run.
% <figuredir> is the directory to write figures to.  this directory should end in
%   an actual directory name because we automatically extract the directory name D
%   and use the prefix 'D_' for figure files.
% <n> is the power of two output matrix size (e.g. 64).  must be greater than 1.
% <s> is the number of slices.  must be greater than 1.
% <llbands> is 1 x 4 with the bandwidths to use in the local linear regression.
%   the four values indicate the bandwidth of the kernel in matrix units for 
%   each of the four dimensions [<n> <n> <s> length(files)] (see localregression4d.m).
%   for example, a value of 2 indicates that the kernel should fall to 0 when we
%   are 2 matrix units away from the center point.  you probably want to use 
%   isotropic kernels, so you should specify different bandwidths for each of the 
%   first three dimensions if voxels are anisotropic.  note that all bandwidths 
%   should be greater than 1 (this is because we need more than one data 
%   point along each dimension in order to perform local linear regression).  also,
%   note that bandwidths can be Inf (this, in effect, results in pure linear regression
%   along the corresponding dimension).  if <files> matches only one run, the 
%   fourth element in <llbands> is ignored, and we perform local linear regression only 
%   over the first three dimensions.
% <mode> (optional) is 
%   0 means normal operation
%   1 means omit the local linear regression and simply use the original fieldmaps.
%     in this case, <llbands> and <dfactor> have no effect.
%   2 means average fieldmaps across runs and use the average fieldmap for all reconstructions.
%     (thus, the local linear regression is omitted.)  in this case, <llbands> and
%     <dfactor> have no effect.
%   3 means average fieldmaps across runs and then subject the result to <llbands>.
%   default: 0.
% <numsteps> (optional) is
%   N means the number of discrete B0 values to use in the reconstruction.  must be at least 2.
%   [N X] means try to use N discrete values (must be at least 2) but do not go below X Hz
%     of spacing between successive values.  so, we may end up using less than N discrete values
%     (but note that we always use at least 2).  in this mode, N can be Inf.
%   default: 50.
% <numframes> (optional) is
%   [] means normal operation
%   N where N is a positive integer means to call grecons with "-f N" when
%     performing reconstructions.  the point is to decrease execution time.
%   default: [].
% <todo> (optional) is a vector of file indices to perform actual reconstructions for.
%   indices refer to the files that are matched by <files>.  <todo> can be 0, which
%   means to do none of the reconstructions.  default: 1:S where S is the number of files
%   matched by <files>.
% <dfactor> (optional) is the factor (positive number) to multiply the bandwidth with to 
%   obtain the granularity at which we perform the local linear regression (and then use cubic 
%   interpolation to obtain the values that we actually want).  (we force the granularity to be at
%   least 3 and no more than the original granularity (as long as it is at least 3).)  
%   if 0, do not perform this speed-up and do the local linear regression directly.  
%   for near-perfect results, you can use 1/3.  for faster execution, try 1/2.  default: 0.
% <isinout> (optional) is whether the data are spiral in-out.  default: 0.
% <gain> (optional) is a number for image scaling (passed to grecons).  default: 1.
% <fvol> (optional) is the "first" volume to inspect.  default: 6.
% <dchack> (optional) is {A B} where A is a fraction of the 99th percentile of the mean brain
%   and B is a vector of DC offsets (one per run) to impose on the final fieldmap.
%   default is [] which means do nothing special.
% <tweaks> (optional) is a non-empty vector of integers indicating tweaks to the final fieldmap
%   lookup indices.  note that the final time-series written to the .mag file corresponds to the first
%   entry in <tweaks>.  note that processing many entries in <tweaks> does not involve much
%   execution time, since a single set of reconstructions is used to generate many tweaking levels.
%   also, note that <dchack> and <tweaks> can co-exist peacefully (<dchack> is handled first, 
%   and then <tweaks> takes effect after that).  can also be a cell vector with as many elements as
%   there are functional runs (each element should be a non-empty vector of integers).  this case
%   allows different sets of tweak values for different runs.
%   default: [0].
% <usesgerun> (optional) is whether to attempt to use sgerun.m to parallel-process the
%   reconstruction of different runs.  if 1, SGE jobs will be submitted and this function
%   will immediately return.  default: 0.
%
% (1) we use grecons to determine the three parameters foff, wx, and wy, 
% which are defined for each slice and each functional run.  (we obtain these 
% values by grep and regexp on the verbose output of grecons.)  summary figures
% are written to parameter1.png, parameter2.png, and parameter3.png (corresponding to 
% foff, wx, and wy).
%
% (2) we then use grecons to calculate the fieldmaps and the brains associated 
% with these fieldmaps.  we write these out for 
% inspection purposes to fieldmap_run%03d.png and brain_run%03d.png.  the same scaling 
% is used for all fieldmaps: [-X,X] where X is the 99th percentile of the absolute 
% value of all pixels in all fieldmaps.  the same scaling is used for all brains: 
% [0,X] where X is the 99th percentile of all pixels in all brains.  note that these
% scaling ranges are also used for all other inspection figures.  we also write out
% a histogram figure for each fieldmap to fieldmap_run%03d_hist.png.
%
% (3) we then calculate the mean across fieldmaps and the mean across fieldmap brains.
% figures are written to fieldmapmean.png and brainmean.png.
%
% (4) assuming that <mode> is 0 or 3, we perform local linear regression to get
% smooth fieldmaps (in both space and time).  we use the values in the 
% fieldmap brains as weights in the local linear regression (note that the
% values are conveniently non-negative).  voxels for which local linear
% regression fails to find a value (i.e. the value gets returned as NaN) are
% explicitly set to 0.  (if <dfactor> is not 0, we write out the intermediate
% coarse local linear approximations to coarse_run001.png, etc.)
%
% (4b) just before finalize the fieldmaps, we apply the <dchack> if necessary.  what this 
% does is to adjust the DC of each fieldmap: considering only those voxels that are greater
% than a certain fraction of the 99th percentile of the mean brain, we enforce a certain
% DC offset (i.e. mean across voxels) on each run's fieldmap.
%
% (5) we then write the fieldmaps to be used to final_run001.png, etc.  then we calculate a 
% robust range of the fieldmaps (see robustrange.m), and discretize the fieldmaps 
% according to <numsteps>.  (we clip voxels whose values fall beyond this range.)  
% we summarize this discretization to discretization.png.  also, we write out the 
% discretized and clipped fieldmap to finalbin_001.png, etc., for inspection purposes.
%
% (6) finally, we reconstruct the functional runs specified by <todo>.  if <numframes>
% is specified, we process only a few volumes.  the reconstruction proceeds by a 
% multi-frequency strategy wherein we calculate many different reconstructions,
% one reconstruction for each discretized B0 value determined in step 5.  for each
% voxel, we determine its value V in the final discretized fieldmap and choose for 
% that voxel's value the value present in the reconstruction that assumed a B0 value
% of V.  the time-series data for all functional runs are written to .mag files 
% (to the same directories that contain the original .7 files).  we also write 
% out the <fvol>th and last volumes for inspection purposes to volumes/run001_tw001a.png, 
% volumes/run001_tw001b.png, and so forth.  the first number indicates the run number, and
% the second number indicates the tweak number.  (if there are less than <fvol> volumes in a 
% functional run, we write out the first and last volumes of that run.)  finally, for the first 
% functional run only, we write out the different B0 reconstructions of the last volume
% to a directory named 'b0' (this is useful for seeing intuitively what effect the different B0 
% values have).
%
% (7) before exiting, we save variables (excluding a few big ones) to files named record*.mat.
% these files are for development purposes only.
%
% please see the code for the exact flags used in the grecons utility.
%
% note that grecons creates files in the directory it is called from as well as in
% the directories that contain the .7 files.  be sure that you do not have two copies 
% of grecons that are trying to reconstruct the same .7 files.  also, we make no attempt
% to clean up the various intermediary files.
%
% some notes on how the case of spiral-in-out data (<isinout>) is handled:
% - directory 'b0' will have two sets of B0 reconstructions, first the spiral-in 
%   volumes and then the spiral-out volumes.
% - the .mag files will have two sets of volumes, first the spiral-in volumes
%   and then the spiral-out volumes (if <numframes> is specified, there will be 
%   2*<numframes> volumes.)
% - inspection volumes will be written to ...a.png, ...b.png, ...c.png, ...d.png,
%   where the first two are the <fvol>th and last volumes of the spiral-in case
%   and the last two are the <fvol>th and last volumes of the spiral-out case.
% - the fieldmap will reflect whichever version (spiral-in or spiral-out) that 
%   grecons outputs.  i think grecons computes the field map from the spiral-out case.
%
% there are several different approaches that you can take towards reconstruction.
% we summarize here several useful approaches:
%   greconsmultifreq.m:
%     1. smooth fieldmaps in space and time.
%     2. average fieldmaps across runs and then smooth in space.
%     3. process each run separately (this is like omitting smoothing over time).
%     4. for a given run, smooth fieldmap in space and then force fieldmap to have 
%        a specific B0 offset (using <dchack>).
%     5. for any of the approaches 1-4, apply a final global tweak to the finalmap 
%        lookup indices (using <tweaks>).  (this provides similar capabilities as 
%        <dchack>, but does so in a different way.)
%   greconsevaluateb0.m:
%     6. evaluate a range of B0 values (this is useful for determining a good 
%        value for <dchack>; see 4 above).
%
% see also greconsevaluateb0.m.
%
% the local regression uses parfor, so to speed up execution you could
% issue matlabpool before calling greconsultifreq.m.
%
% random sysadmin notes:
% - you may have to do "yum install glibc.i686" in order to be able to execute grecons.
%
% history:
% 2010/10/11 - remove <externalfmap>; implement usesgerun
% 2010/10/01 - implement <externalfmap>
% 2010/07/09 - nicer color order for parameter plots
% 2010/06/24 - allow <tweaks> to be a cell vector
% 2010/06/10 - tweak some names for figure files
% 2010/05/29 - now, write out histogram figure for each individual raw fieldmap.
%              new special mode for <numsteps>.
% 2010/05/27 - implement <tweaks>.
% 2010/05/26 - initial version (no more thresh or brain mask; require grecons rev 95; 
%              faster local regression; write out coarse local regression results;
%              implement spiral in-out case).  2:15 pm.

% TODO: would line plots of fieldmaps be useful?
% TODO: should we implement local regression degree 0?
% TODO: make this truly parallelized
% **TODO: write out the motion corrected volumes for free!!!  auto correct.  coregistervolumes.m
% TODO: automatic tweaking? [lame]
% **TODO: figure with mean over slice  [this very important for figuring out temporal bandwidth].
%         with some smoothing localregression: hold on; plot(1:12,localregression(1:12,xx,1:12,1,[],2.5),'mo-')
%         >> good = mean(a.brains,4) > prctile(a.brains(:),99)*.1;
%         >> good = find(good);
%         >> check=mean(subscript(squish(a.fieldmaps,3),{good ':'}),1);

% internal constants
dformat = 'int16';                       % the data format used by grecons
bname = 'brain_run%03d';                 % file names for the brains
bname0 = 'brainmean';                    % file names for the brain mean
fname = 'fieldmap_run%03d';              % file names for the fieldmaps
fnameh = 'fieldmaphist_run%03d';        % file name for the fieldmap histogram
fname0 = 'fieldmapmean';                 % file names for the fieldmap mean
gname0 = 'coarse_run%03d';             % file name for the coarse local linear
zname1 = 'final_run%03d';                % file names for the final
zname2 = 'finalbin_run%03d';             % file names for the final (binned)
movdir = 'b0';                           % directory name for the movie images
movname = '%06d';                        % file name for the movie images
inamedir = 'volumes';
iname = 'run%03d_tw%03d';                % file names for the inspections
pnames = {'foff' 'wx' 'wy'};             % names of the various parameters
nparams = length(pnames);                % total number of parameters
ptile = 99;                              % percentile for scaling purposes
rname = 'record';                        % the .mat file we write out
dname = 'discretization';                % file name for the discretization figure
parname = 'parameter%d';                 % file name for the parameter figures

% inputs
if ~exist('mode','var') || isempty(mode)
  mode = 0;
end
if ~exist('numsteps','var') || isempty(numsteps)
  numsteps = 50;
end
if ~exist('numframes','var') || isempty(numframes)
  numframes = [];
end
if ~exist('todo','var') || isempty(todo)
  todo = [];  % deal with later
end
if ~exist('dfactor','var') || isempty(dfactor)
  dfactor = 0;
end
if ~exist('isinout','var') || isempty(isinout)
  isinout = 0;
end
if ~exist('gain','var') || isempty(gain)
  gain = 1;
end
if ~exist('fvol','var') || isempty(fvol)
  fvol = 6;
end
if ~exist('dchack','var') || isempty(dchack)
  dchack = [];
end
if ~exist('tweaks','var') || isempty(tweaks)
  tweaks = [0];
end
% if ~exist('externalfmap','var') || isempty(externalfmap)
%   externalfmap = [];
% end
if ~exist('usesgerun','var') || isempty(usesgerun)
  usesgerun = 0;
end

% calc
prefix = [stripfile(figuredir,1) '_'];
nv = choose(isinout,4,2);  % how many inspection volumes will we have at the end
letters = makeletters(nv);  % e.g. {'a' 'b' 'c' 'd'}

% make figuredir if necessary
mkdirquiet(figuredir);

% match the filenames
files = matchfiles(files);
assert(~isempty(files),'<files> does not match at least one file');
% if ~isempty(externalfmap)
%   externalfmap = matchfiles(externalfmap);
%   assert(length(externalfmap)==length(files),'number of <files> and number of <externalfmap> do not match');
% end
fprintf('we are processing these files:\n'); fprintf(' %s\n',files{:});

% more input handling
if ~iscell(tweaks)
  tweaks = repmat({tweaks},[1 length(files)]);
end
assert(length(tweaks)==length(files));

% get fieldmap parameters
params = zeros(s,nparams,length(files));
for p=1:length(files)
  for q=1:s
    result = unix_wrapper(sprintf('%s -s %d -e %d -v -B 2 -i 1 -G %.10f -N %d -O %s | grep fof',grecons,q,q,gain,n,files{p}));
    tokens = regexp(result,'.+=(?<foff>.+?) Hz.+?wx, wy =(?<wx>.+?) (?<wy>.+?) rad','names');
    params(q,:,p) = [str2num(tokens.foff) str2num(tokens.wx) str2num(tokens.wy)];
  end
end

% inspect parameters
for p=1:nparams
  figureprep; hold all;
  set(gca,'ColorOrder',jet(length(files)));
  h = plot(squish(params(:,p,:),2));
  legend(h,cellfun(@(x) stripfile(x,1),files,'UniformOutput',0));
  xlabel('slice number'); ylabel('parameter value');
  title(sprintf('parameter %s',pnames{p}));
  figurewrite(sprintf([prefix parname],p),[],[],figuredir);
end

% load in fieldmaps (PHASE) and brains (MAG)
fieldmaps = zeros(n,n,s,length(files));
brains = zeros(n,n,s,length(files));
% if isempty(externalfmap)
%   fmapfiles = files;
% else
%   fmapfiles = externalfmap;
% end
fmapfiles = files;
for p=1:length(fmapfiles)
  for q=1:s
    fieldmaps(:,:,q,p) = loadbinary(sprintf([fmapfiles{p} '.B0.%03d'],q),dformat,[n n]);
    brains(:,:,q,p) = loadbinary(sprintf([fmapfiles{p} '.%03d'],q),dformat,[n n]);
  end
end
    % apparently, the grecons earlier already made the files, so this isn't necessary:
    %  assert(unix(sprintf('%s -B 2 -i 1 -N %d -O %s',grecons,n,files{p}))==0);
% for rev100 !!!
%   for p=1:length(externalfmap)
%     temp = loadbinary(externalfmap{p},dformat,[n n 2*s]);
%     brains(:,:,:,p) = temp(:,:,1:2:end);
%     fieldmaps(:,:,:,p) = temp(:,:,2:2:end);
%   end
%end

% write fieldmaps and brains to .png files for inspection
fmx = prctile(abs(fieldmaps(:)),ptile);
bmx = prctile(brains(:),ptile);
for p=1:length(files)
  imwrite(uint8(255*makeimagestack(fieldmaps(:,:,:,p),[-fmx fmx])),jet(256),sprintf([figuredir '/' prefix fname '.png'],p));
  imwrite(uint8(255*makeimagestack(brains(:,:,:,p),[0 bmx])),gray(256),sprintf([figuredir '/' prefix bname '.png'],p));
    % histogram figure
  ftemp = flatten(fieldmaps(:,:,:,p));
  [d,mn,mx] = robustrange(ftemp);
  figureprep; hold on;
  [nn,xx] = hist(ftemp,linspace(mn-(mx-mn)*0.5,mx+(mx-mn)*0.5,100));
  bar(xx,nn,1);
  ax = axis; axis([mn-(mx-mn)*0.6 mx+(mx-mn)*0.6 ax(3:4)]);
  xlabel('fieldmap values');
  figurewrite(sprintf([prefix fnameh],p),[],[],figuredir);
end

% calc mean across runs
fieldmap0 = mean(fieldmaps,4);
brain0 = mean(brains,4);

% write some figures
imwrite(uint8(255*makeimagestack(fieldmap0,[-fmx fmx])),jet(256),[figuredir '/' prefix fname0 '.png']);
imwrite(uint8(255*makeimagestack(brain0,[0 bmx])),gray(256),[figuredir '/' prefix bname0 '.png']);

% prep
[xx,yy,zz,tt] = ndgrid(1:n,1:n,1:s,1:length(files));  % indices of all voxels

% handle special cases
switch mode

% in this case, just use the original fieldmaps for reconstruction
case 1
  final = fieldmaps;

% in this case, just use the mean fieldmap for reconstruction
case 2
  final = repmat(fieldmap0,[1 1 1 length(files)]);

% in this case, calculate the mean fieldmap as a starting point and then allow llbands to take effect
case 3  
  fieldmaps = repmat(fieldmap0,[1 1 1 length(files)]);  % NOTE THAT THIS CLOBBERS fieldmaps !
  brains = repmat(brain0,[1 1 1 length(files)]);        % NOTE THAT THIS CLOBBERS brains !

end

% some checks
assert(all(~isnan(fieldmaps(:))));
assert(all(~isnan(brains(:))));

% perform local linear regression to smooth the fieldmap
if ismember(mode,[0 3])

  % perform regression for each bandwidth and write inspection figures
  if dfactor==0
    if length(files) > 1
      final = localregression4d(xx,yy,zz,tt,fieldmaps,xx,yy,zz,tt,1,'epan',llbands,brains,1);
    else
      final = localregression3d(xx,yy,zz,   fieldmaps,xx,yy,zz,   1,'epan',llbands(1:3),brains,1);
    end
    final(isnan(final)) = 0;
  else
    isize = ceil(([n n s length(files)]-1) ./ (llbands*dfactor));  % number of elements along each dimension for the quick regression
    isize = max([3 3 3 3; min([n n s length(files); max([3 3 3 3; isize],[],1)],[],1)],[],1);  % at least 3, no more than original granularity, at least 3
    if length(files) > 1
      [xx2,yy2,zz2,tt2] = ndgrid(linspace(1,n,isize(1)),linspace(1,n,isize(2)),linspace(1,s,isize(3)),linspace(1,length(files),isize(4)));
      coarse = localregression4d(xx,yy,zz,tt,fieldmaps,xx2,yy2,zz2,tt2,1,'epan',llbands,brains,1);
      coarse(isnan(coarse)) = 0;  % don't let NaNs propagate!
      coarse = clipoutliers(coarse);
      final = [];
      for q=1:length(files)  % to reduce memory usage, let's interpolate each time-slice separately
        final(:,:,:,q) = interpn(linspace(1,n,isize(1)),linspace(1,n,isize(2)),linspace(1,s,isize(3)),linspace(1,length(files),isize(4)), ...
                                  coarse,xx(:,:,:,q),yy(:,:,:,q),zz(:,:,:,q),tt(:,:,:,q),'*cubic');
      end
    else
      [xx2,yy2,zz2] = ndgrid(linspace(1,n,isize(1)),linspace(1,n,isize(2)),linspace(1,s,isize(3)));
      coarse = localregression3d(xx,yy,zz,fieldmaps,xx2,yy2,zz2,1,'epan',llbands(1:3),brains,1);
      coarse(isnan(coarse)) = 0;  % don't let NaNs propagate!
      coarse = clipoutliers(coarse);
      final = interpn(linspace(1,n,isize(1)),linspace(1,n,isize(2)),linspace(1,s,isize(3)),coarse,xx,yy,zz,'*cubic');
    end
    assert(all(~isnan(final(:))));  % HMM, IF GRID IS THE SAME AS ORIGINAL, WE DO UNNECESSARY COMPUTATION..
    for q=1:size(coarse,4)
      imwrite(uint8(255*makeimagestack(coarse(:,:,:,q),[-fmx fmx])),jet(256),sprintf([figuredir '/' prefix gname0 '.png'],q));
    end
  end

end

% deal with dchack
if ~isempty(dchack)
  ok = brain0 > prctile(brain0(:),99)*dchack{1};
  for p=1:size(final,4)
    final(:,:,:,p) = final(:,:,:,p) - mean(subscript(final(:,:,:,p),ok)) + dchack{2}(p);
  end
end

% write out the final field map we are going to use
for p=1:size(final,4)
  imwrite(uint8(255*makeimagestack(final(:,:,:,p),[-fmx fmx])),jet(256),sprintf([figuredir '/' prefix zname1 '.png'],p));
end

% figure out discretization values for the fieldmap
[d,mn,mx] = robustrange(final(:));
if length(numsteps) > 1
  maxnum = floor((mx-mn)/numsteps(2)) + 1;  % this is the maximum number that will not go below numsteps(2) spacing
  numsteps = max(2,min(maxnum,numsteps(1)));  % if numsteps(1) is too large, bring it down to maxnum; but make sure at least 2
  % notice that the above line clobbers numsteps!
end
values = linspace(mn,mx,numsteps);

% summary figure
figureprep; hold on;
[nn,xx] = hist(final(:),linspace(mn-(mx-mn)*0.5,mx+(mx-mn)*0.5,100));
bar(xx,nn,1);
straightline(values,'v','r-');
ax = axis; axis([mn-(mx-mn)*0.6 mx+(mx-mn)*0.6 ax(3:4)]);
xlabel('fieldmap values');
title(sprintf('chopped %.1f percent; number of steps is %d; step is %.1f', ...
  sum(final(:) < mn | final(:) > mx) / length(final(:)) * 100, ...
  numsteps,values(2)-values(1)));
figurewrite([prefix dname],[],[],figuredir);

% calculate indices of how we approximate the field map (integers from 1 to numsteps)
indices = round((final-mn) / (mx-mn) * (numsteps-1)) + 1;
indices(indices<1) = 1;
indices(indices>numsteps) = numsteps;

% write out the approximation to the field map
for p=1:size(indices,4)
  imwrite(uint8(255*makeimagestack((indices(:,:,:,p)-1) / (numsteps-1) * (mx-mn) + mn,[-fmx fmx])),jet(256),sprintf([figuredir '/' prefix zname2 '.png'],p));
end

% save
clear xx yy zz tt xx2 yy2 zz2 tt2;
save([figuredir '/' prefix rname '.mat']);

%%%%%%%%%%% RECONSTRUCTION STAGE

% do the reconstruction
if ~isequal(todo,0)
  todo = choose(isempty(todo),1:length(files),todo);

  % HERE, WE WOULD INSERT THE FOR-LOOP.
  % INSTEAD, LET'S DISPATCH USING SGERUN IF THE USER WANTS IT:

  cmd = loadmulti(strrep(which('greconsmultifreq'),'greconsmultifreq.m','greconsmultifreq.mat'),'cmd');
  if usesgerun
    sgerun(cmd,0,1,1:length(todo));
  else
    for jobindex=1:length(todo)
      eval(cmd);
    end
  end

end

return;

%%%%%%%%%%% FOR-LOOP [FOR DEVELOPMENT ONLY]

cmd = inputmulti

% for p=1:length(todo)
    p = jobindex;

    % calc
    ii = todo(p);  % the actual file index
  
    % reconstruct each slice at each b0 step and collect them up [memory is a problem, so do it slice by slice!]
    inspections = [];
    b0frames = zeros(0,0,dformat);
    finalvol = zeros(0,0,dformat);
    for q=1:s
  
      % suck the data in
      vol = zeros(0,0,dformat);
      for r=1:numsteps
        unix_wrapper(sprintf('%s %s -s %d -e %d -B 0 -o %.10f -i 1 -G %.10f -N %d -O %s', ...
          grecons,choose(isempty(numframes),'',sprintf('-f %d',numframes)),q,q, ...
          values(r),gain,n,files{ii}));
        temp = loadbinary([files{ii} '.mag'],dformat,[n n 0]);
        if isinout && ~isempty(numframes)
          vol(:,:,:,1,r) = temp(:,:,[1:numframes end-numframes+1:end]);  % got to ignore the middle zeros
        else
          vol(:,:,:,1,r) = temp;  % X x Y x time x 1 x steps
        end
        % ok, when <isinout>, vol will have two sets of volumes
      end
      
      % record things we care about
      if p==1
        if isinout
          b0frames(:,:,q,1:numsteps) = reshape(vol(:,:,end/2,1,:),n,n,1,[]);  % X x Y x Z x steps
          b0frames(:,:,q,numsteps+(1:numsteps)) = reshape(vol(:,:,end,1,:),n,n,1,[]);  % X x Y x Z x steps
        else
          b0frames(:,:,q,:) = reshape(vol(:,:,end,1,:),n,n,1,[]);  % X x Y x Z x steps
        end
      end
      vol = squish(permute(vol,[5 1 2 4 3]),4);  % 100*X*Y*1 x time
      for r=1:length(tweaks{ii})
        vol2 = vol(flatten(max(1,min(numsteps,indices(:,:,q,ii) + tweaks{ii}(r)))) + numsteps * (0:n*n-1),:);  % X*Y*1 x time    [TODO: SHOULD THIS BE A FUNCTION?]
        finalvol(:,:,:,q,r) = permute(reshape(vol2,n,n,1,[]),[1 2 4 3]);  % X x Y x time x Z x tweaks
      end
      
    end
    
    % show the last volume at all the b0 values (only the first case though)
    if p==1
      mkdirquiet([figuredir '/' prefix movdir]);
      for q=1:size(b0frames,4)  % use size of b0frames because it may be double the number of steps
        imwrite(uint8(255*makeimagestack(double(b0frames(:,:,:,q)),[0 bmx])),gray(256),sprintf([figuredir '/' prefix movdir '/' movname '.png'],q));
      end
    end
    
    % save data to the binary file
    savebinary([files{ii} '.mag'],dformat,finalvol(:,:,:,:,1));  % only the first tweak gets saved
    
    % write inspection volumes
    mkdirquiet([figuredir '/' prefix inamedir]);
    for r=1:length(tweaks{ii})
      if isinout
        fvol0 = choose(size(finalvol,3)/2 < fvol,1,fvol);
        inspections = cat(4,inspections,permute(double(finalvol(:,:,[fvol0 end/2 end/2+fvol0 end],:,r)),[1 2 4 3 5]));
      else
        fvol0 = choose(size(finalvol,3) < fvol,1,fvol);
        inspections = cat(4,inspections,permute(double(finalvol(:,:,[fvol0 end],:,r)),[1 2 4 3 5]));
      end
      for q=1:nv
        imwrite(uint8(255*makeimagestack(inspections(:,:,:,end-q+1),[0 bmx])),gray(256),sprintf([figuredir '/' prefix inamedir '/' iname letters{end-q+1} '.png'],p,r));
      end
    end

    % save record
    clear b0frames finalvol vol vol2;
    save(sprintf([figuredir '/' prefix rname '%03d.mat'],p));
  
%  end

save('greconsmultifreq.mat','cmd');





% This is the documentation returned by grecons14 rev 85:
%
% Usage: grecons14 [-s isl1] [-e isl2] [-f nframes] [-B b0map] [-o freq_offset]
%        [-i imgtyp] [-N npix] [-g gamma] [-G gain] [-z zoomf] [-F fermif]
%        [-n nav] [-bcmOuv] [-C coil_list] [-x wxoff] [-y wyoff] raw_file
%        isl1, isl2 are start and end slices to do
%        b0map = 0 for no shim correction
%              = 1 for normal correction [default]
%              = 2 for Bo map only
%        freq_offset (Hz) adds to b0map value
%        imgtyp is a bit map: 1=mag, 2=phase, 4=real, 8=imag
%        npix is power of two output matrix size
%        gamma is in Hz/G [default 4257] 
%        gain is factor for image scaling [default 1] 
%        zoomf is magnification factor
%        fermif is filter radius/kmax, 0 for none [0]
%        nav = 1 for navigator correction [(views==1)? 1:0]
%        -b to swap bytes in output images
%        -c to save coil images in multicoil
%        -m for no maxwell field correction
%        -O to generate images in single output file
%        -u to generate images for unfold
%        -v to print slice/echo notifications (verbose is slower)
%        coil_list is concatenated list of coils to use,
%            e.g., -C 0147 for coils 0, 1, 4 & 7
%        wxoff is offset for field map xslope
%        wyoff is offset for field map yslope
% This is source rev 85: 10/5/2008  signa rev 14.5


% This is the documentation returned by grecons14 rev 95:
%
% Usage: grecons14_rev95 [-s isl1] [-e isl2] [-f nframes] [-B b0map] [-o freq_offset]
%        [-i imgtyp] [-N npix] [-g gamma] [-G gain] [-z zoomf] [-F fermif]
%        -bcmnOuv] [-C coil_list] [-x wxoff] [-y wyoff] raw_file
%        isl1, isl2 are start and end slices to do
%        b0map = 0 for no shim correction
%              = 1 for normal correction [default]
%              = 2 for Bo map only
%        freq_offset (Hz) adds to b0map value
%        imgtyp is a bit map: 1=mag, 2=phase, 4=real, 8=imag
%        npix is power of two output matrix size
%        gamma is in Hz/G [default 4257] 
%        gain is factor for image scaling [default 1] 
%        zoomf is magnification factor
%        fermif is filter radius/kmax, 0 for none [0]
%        -b to swap bytes in output images
%        -c to save coil images in multicoil
%        -m for no maxwell field correction
%        -n for no navigator correction [(views==1)? 1:0]
%        -O to generate images in single output file
%        -u to generate images for unfold
%        -v to print slice/echo notifications (verbose is slower)
%        coil_list is concatenated list of coils to use,
%            e.g., -C 0147 for coils 0, 1, 4 & 7
%        wxoff is offset for field map xslope
%        wyoff is offset for field map yslope
% This is source rev 95: 5/22/10  signa rev 14.5





%%%%%%%%%% INTERNAL NOTES

% <externalfmap> (optional) is a file pattern (see matchfiles.m) matching one or more
%   .7 files.  [[####each file should be int16 and have dimensions <n> x <n> x 2*<s>, where the
%   odd slices are the magnitude brains and the even slices are the fieldmaps.
%   the number of .B0 files should match the number of <files>.]]  if supplied, we use
%   these data for the fieldmaps instead of what is contained in the functional runs.
%   default is [] which means do nothing special.

% (2) we then use grecons to calculate the fieldmaps and the brains associated 
% with these fieldmaps (unless <externalfmap> is specified, in which case we load
% the fieldmaps and brains from the specified files).

function [subj] = AG_statmap_random(subj,data_patname,regsname,selname,new_map_patname,extra_arg)
		  
% Creates statmap at randomn
%
% [SUBJ] = AG_statmap_random(SUBJ,DATA_PATNAME,REGSNAME,SELNAME,NEW_MAP_PATNAME,EXTRA_ARG)
%
%
% IMPORTANT: The regressors matrix must be a single row vector in
% order for such calculations to occur.
%
% Adds the following objects:
% - statmap pattern object
%
% Only uses those TRs labelled with a 1 in the SELNAME selector.
%
% EXTRA_ARG is ignored.
%
% License:
%=====================================================================
%
% This is part of the Princeton MVPA toolbox, released under
% the GPL. See http://www.csbmb.princeton.edu/mvpa for more
% information.
% 
% The Princeton MVPA toolbox is available free and
% unsupported to those who might find it useful. We do not
% take any responsibility whatsoever for any problems that
% you have related to the use of the MVPA toolbox.
%
% ======================================================================

defaults.extra_arg = [];
defaults.cur_iteration = NaN;

% default arguments
args = propval({extra_arg}, defaults); 

pat  = get_mat(subj,'pattern',data_patname);
regs = get_mat(subj,'regressors',regsname);
sel  = get_mat(subj,'selector',selname);

%sanity_check(pat,regs,sel,args);

TRs_to_use = find(sel==1);

% Note: don't forget to exclude rest timepoints, unless your
% function definitely requires them

pat   = pat(:,TRs_to_use);
regs = regs(:,TRs_to_use);

if isempty(pat)
  error('Cannot compute statmap on an empty pattern');
end
  
for i = 1:size(pat,1)
    if  isequal(diff(pat(i,:)), zeros(1,size(pat,2)-1))
        constant_rows(i)=1;
    else
        constant_rows(i)=0;
    end
end

%generate a random pattern the size of 'pat'
new_pat = randn(size(pat,1),1);

% Create the new pattern, using the same mask as before
masked_by = get_objfield(subj,'pattern',data_patname,'masked_by');
subj = initset_object(subj, 'pattern', new_map_patname, new_pat, 'masked_by', ...
		      masked_by);

% add bookkeeping history
hist = sprintf('Created by %s',mfilename());
subj = add_history(subj,'pattern',new_map_patname,hist);

created.function = mfilename();
created.data_patname = data_patname;
created.regsname = regsname;
created.selname = selname;
created.new_map_patname = new_map_patname;
created.args = args;
subj = add_created(subj,'pattern',new_map_patname,created);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xcorr] = statmap_xcorr_logic(pat,regs,args)

% check for C language version
if exist('compute_xcorr')
  xcorr = compute_xcorr(regs', pat');
else
  
  warning(['compute_xcorr.c has not been compiled. Computation will ' ...
           'be slow.']);
  
  nVoxels = size(pat, 1);

  xcorr = zeros(nVoxels,1);

  for v = 1:nVoxels   

    %compute correlation coefficient between voxel, regressor
    if (std(regs) == 0 || std(pat(v,:)) == 0)
      error(['Cannot compute cross-correlation if std of input is zero. ' ...
             'Check your input for constants.']);
    end  
    xcorr(v) = 1 - pdist([regs; pat(v,:)], 'correlation');
  end 

end

% take the absolute value: negatively correlated voxels are still informative!
xcorr = abs(xcorr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = sanity_check(pat,regs,sel,args)

if size(regs,1)>1
  error(sprintf('Regressor is not a row vector.', regs));
end

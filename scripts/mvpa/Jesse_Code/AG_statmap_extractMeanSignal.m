function meanSignal = AG_statmap_extractMeanSignal(subj,data_patname,regsname,selname,new_map_patname,extra_arg)

% Use the anova to select features that vary between conditions
%
% [SUBJ] = STATMAP_ANOVA(SUBJ,DATA_PATNAME,REGSNAME,NEW_MAP_PATNAME,EXTRA_ARG);
%
% Adds the following objects:
% - statmap pattern object
%
% Updates the subject structure by creating a new pattern,
% NEW_MAP_PATNAME, that contains a vector of P-values from the ANOVA.
%
% Uses all the conditions in REGSNAME. If you only want to use a
% subset of them, create a new regressors object with only those
% conditions
%
% Only uses those TRs labelled with a 1 in the SELNAME selector,
% and where there's an active condition in the REGSNAME regressors matrix
%
% There should be functionality in here for return the F values as
% well (e.g. an optional argument 'MAP_TYPE' that defaults to 'p') xxx
%
% All statmap functions have to take in an EXTRA_ARG argument from
% FEATURE_SELECT.M. In this case, it has only one optional field:
%
% - USE_MVPA_VER (optional, default = false). If true, this will use
% the MVPA anova function (ANOVA1_MVPA.M) rather than the Stats
% toolbox ANOVA1.M

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



pat  = get_mat(subj,'pattern',data_patname);
regs = get_mat(subj,'regressors',regsname);
sel  = get_mat(subj,'selector',selname);


TRs_to_use = find(sel==2);

pat   = pat(:,TRs_to_use);
regs = regs(:,TRs_to_use);


meanSignal = mean(pat);






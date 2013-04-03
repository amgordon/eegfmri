function PM_GroupWholeshebang(subpar, flags)

if isstruct(subpar) % if it is par_params struct
    par = subpar;
else % assume subject string
    par = EF_GroupParams(subpar);
end

if ismember('m', flags);
            EF_groupModelSpec(par);
end

if ismember('e', flags);
    for tsk = 1:length(par.tasks)
        for cnd = 1:length(par.task{tsk}.cons)
            
            EF_group_mod_est(par, tsk, cnd);
            
        end
    end
end

if ismember('c', flags);
    for tsk = 1:length(par.tasks)
        for cnd = 1:length(par.task{tsk}.cons)
            EF_groupsetcontrasts(par, tsk, cnd);            
        end
    end
end
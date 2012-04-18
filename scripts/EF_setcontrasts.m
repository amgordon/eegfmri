function PM_setcontrasts(subpar)
% parSetContrasts(subpar)
% Sets contrasts.
% 
% subpar refers to either the subject number (in string notation) or
% the par struct generated by par = par_params(subject)
% 
% made from av_setContrasts  1/28/08 jbh
% modified to subpar format, etc, 7/24/08 jbh



origdir = pwd;


% ---load par params if need be---
if isstruct(subpar) % if it is par_params struct
    par = subpar;
else % assume subject string
    par = par_params(subpar);
end

STAT = par.constat;

cd(par.analysisdir);
fprintf('\nLoading SPM...');
load SPM
fprintf('Done');

Xsize = size(SPM.xX.xKXs.X,2);

padCon = @padConWithZeros;

regNames = SPM.xX.name;

% T-contrasts
%---------------------------------------------------------------------------



% so you want to make a contrast...
% instead of manually typing out the contrast matrix, I'll create them
% based on contrast names.  Anything 'fancier' than simple contrasts (e.g.
% parametric contrasts) will be specified after this section.  This will
% also produce the inverse contrasts of anything specified.  Note that for
% balancing ease of interpretation with shorthand, all conditions must
% contain only one initial cap [A-Z]
% current conditions:


idx.hits = ~cellfun('isempty', strfind(regNames,'hits'));
idx.CRs = ~cellfun('isempty', strfind(regNames,'CRs'));
idx.HC = ~cellfun('isempty', strfind(regNames,'HC'));
idx.ERP = ~cellfun('isempty', strfind(regNames,'ERPAmp'));
idx.all = ~cellfun('isempty', strfind(regNames,'all'));

idx.finger1 = ~cellfun('isempty', strfind(regNames,'finger1'));
idx.finger2 = ~cellfun('isempty', strfind(regNames,'finger2'));
idx.finger3 = ~cellfun('isempty', strfind(regNames,'finger3'));
idx.finger4 = ~cellfun('isempty', strfind(regNames,'finger4'));
idx.finger5 = ~cellfun('isempty', strfind(regNames,'finger5'));
idx.finger = ~cellfun('isempty', strfind(regNames,'finger'));

switch par.ONAnalysisType
    case 'hitsVsCRs'
        con.hitsVsCRs = idx.hits - idx.CRs;
    case 'hitsVsCRs_HC'
        con.hitsVsCRs_HC = idx.HC .* (idx.hits - idx.CRs);
    case 'hitsVsCRs_byERP'
        con.hitsByERP = idx.hits .* idx.ERP;
        con.CRsByERP = idx.CRs .* idx.ERP;
        con.hitsAndCRsByERP = (idx.hits + idx.CRs) .* idx.ERP;
        con.hitsVsCRs_X_ERPAmp = (idx.hits - idx.CRs) .* idx.ERP;
    case 'hitsVsCRs_byERP_HC'
        con.hitsByERP_HC = idx.hits .* idx.ERP .* idx.HC;
        con.CRsByERP_HC = idx.CRs .* idx.ERP .* idx.HC;
        con.hitsAndCRsByERP_HC = (idx.hits + idx.CRs) .* idx.ERP .* idx.HC;
        con.hitsVsCRs_X_ERPAmp_HC = (idx.hits - idx.CRs) .* idx.ERP .* idx.HC;
    case 'alltrials_byERP'
        con.allByERP = idx.all .* idx.ERP;
    case 'alltrials_byERP_HC'
        con.allByERP = idx.all .* idx.ERP .* idx.HC;
    case 'buttonPresses_byERP_RTLocked'
        con.allByERP = idx.ERP;
        con.finger1ERP = idx.ERP .* idx.finger1;
        con.finger2ERP = idx.ERP .* idx.finger2;
        con.finger3ERP = idx.ERP .* idx.finger3;
        con.finger4ERP = idx.ERP .* idx.finger4;
        con.finger5ERP = idx.ERP .* idx.finger5;
        
        con.all_events_vs_fix = idx.finger .* ~idx.ERP;
        con.finger1 = idx.finger1 .* ~idx.ERP;
        con.finger2 = idx.finger2 .* ~idx.ERP;
        con.finger3 = idx.finger3 .* ~idx.ERP;
        con.finger4 = idx.finger4 .* ~idx.ERP;
        con.finger5 = idx.finger5 .* ~idx.ERP;
    case 'buttonPresses_byERP'
        con.allByERP = idx.ERP;
    case 'buttonPress_by_spectralPower'
        con.allByERP = idx.ERP;
        con.finger1ERP = idx.ERP .* idx.finger1;
        con.finger2ERP = idx.ERP .* idx.finger2;
        con.finger3ERP = idx.ERP .* idx.finger3;
        con.finger4ERP = idx.ERP .* idx.finger4;
        con.finger5ERP = idx.ERP .* idx.finger5;
        
        con.all_events_vs_fix = idx.finger .* ~idx.ERP;
        con.finger1 = idx.finger1 .* ~idx.ERP;
        con.finger2 = idx.finger2 .* ~idx.ERP;
        con.finger3 = idx.finger3 .* ~idx.ERP;
        con.finger4 = idx.finger4 .* ~idx.ERP;
        con.finger5 = idx.finger5 .* ~idx.ERP;
        
    case 'buttonPress_by_spectralPower_RTLocked'
        
        con.allByERP = idx.ERP;
        con.finger1ERP = idx.ERP .* idx.finger1;
        con.finger2ERP = idx.ERP .* idx.finger2;
        con.finger3ERP = idx.ERP .* idx.finger3;
        con.finger4ERP = idx.ERP .* idx.finger4;
        con.finger5ERP = idx.ERP .* idx.finger5;
        
        con.all_events_vs_fix = idx.finger .* ~idx.ERP;
        con.finger1 = idx.finger1 .* ~idx.ERP;
        con.finger2 = idx.finger2 .* ~idx.ERP;
        con.finger3 = idx.finger3 .* ~idx.ERP;
        con.finger4 = idx.finger4 .* ~idx.ERP;
        con.finger5 = idx.finger5 .* ~idx.ERP;
end


fn_con = fieldnames(con);

for f = 1:length(fn_con)
    cnames{f} = fn_con{f};
    cvals{f} = con.(fn_con{f});
end


% preallocate
con_name(1:length(cnames)) = {''};
con_vals = cell(1, length(cnames));

%con_vals(1:length(cnames)) = {zeros(1,length(cnames))};

for Tt = 1:length(cnames)

    % make names
    con_name{Tt} = cnames{Tt};

    con_vals{Tt} = double(cvals{Tt});

end



% put contrasts into SPM/write to file

fprintf('\nBeginning contrasts on subject %s\n', par.substr);


cfid = fopen('conlist','wt');
fprintf(cfid, 'Contrasts for Sub %s\nLast run on %s\n', par.substr, date);

% Loop over created contrasts
%-------------------------------------------------------------------
for k=1:length(con_vals)

    % Basic checking of contrast
    %-------------------------------------------------------------------
    [c,I,emsg,imsg] = spm_conman('ParseCon',con_vals{k},SPM.xX.xKXs,STAT);
    if ~isempty(emsg)
        disp(emsg);
        error('Error in contrast specification');
    else
        disp(imsg);
    end;

    % Fill-in the contrast structure
    %-------------------------------------------------------------------
    if all(I)
        DxCon = spm_FcUtil('Set',con_name{k},STAT,'c',c,SPM.xX.xKXs);
    else
        DxCon = [];
    end

    % Append to SPM.xCon. SPM will automatically save any contrasts that
    % evaluate successfully.
    %-------------------------------------------------------------------
    if isempty(SPM.xCon)
        SPM.xCon = DxCon;
    elseif ~isempty(DxCon)
        SPM.xCon(end+1) = DxCon;
    end
    SPM = spm_contrasts(SPM,length(SPM.xCon));
        
    fprintf(fopen('conlist','at'),'%d: %s\n%s\n\n',k, con_name{k},num2str(con_vals{k}));
end

fclose(cfid);
copyfile('conlist',[par.logdir filesep 'conlist-' date]);

% Change back directory
cd(origdir);

fclose('all')
return;

% put in here, so I don't have to add it to path or go back to scripts dir
% to execute...
function con = padConWithZeros( cIn, Xsize )

conLength = length(cIn);
nZeros = Xsize - conLength;
con = [cIn zeros(1,nZeros)];


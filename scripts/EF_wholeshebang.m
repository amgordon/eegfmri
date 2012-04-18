
function EF_wholeshebang(subpar, flags)
% function wholeshebang(subpar, flags)
% runs an entire subject start to finish.
%
% subpar refers to either the subject number (in string notation) or
% the par struct generated by par = par_params(subject)
%
% flags defaults to all stages, or individual stages can be chosen from:
% m = makevols
% a = make anatomicals
% s = slice timing
% l = realignment
% c = coregister inplane anat to mean func, and hires to inplane
% g = segment and normalize gray matter
% n = normalize functionals
% h = smooth functionals
% k = make specmask
% z = make artifact indices
% o = make and move onsets
% r = make mvpa regressors
% v = make mvpa regressors
% p = specify model
% e = run model
% t = make contrasts
%
% M = mail .ps file to yourself
% 
% modified from older version by jbh on 7/23/08

origdir = pwd;

if ~exist('subpar', 'var')
    error('Must specify subject!');
end

% ---load par params if need be---
if isstruct(subpar) % if it is par_params struct
    par = subpar;
else % assume subject string
    par = EF_params(subpar);
end

if ~exist('flags', 'var')
    fprintf('\n\n***WARNING: Assuming default flags: doing everything!***\n\n');
    pause(5); % built-in bailout time for if you really didn't want to take the default option
    flags = 'maslcgnhkorpetM';
end

%set up logdir if it doesn't exist (for saving par to, etc)
if ~exist(par.logdir, 'dir'); mkdir(par.logdir); end
%save par file so that you know what you did!
parfile = [par.logdir filesep ['par-' date '.mat']];
save(parfile, 'par');


% Preprocessing
% Assumes 'raw' directory where 4d scans are named 'scan0x'; behavioral
% directory; and anatomical directory
if ismember('i',flags); EF_make3dNiftis(par); end

if ismember('s',flags); PM_slicetime(par); end

if ismember('l',flags); PM_realign(par); end

if ismember('c',flags); PM_coregwrapper(par); end

if ismember('g',flags); PM_segnorm(par); end

if ismember('n',flags); EF_normfuncs(par); end

if ismember('h',flags); PM_smoothfuncs(par); end

if ismember('w',flags); PM_normSPGR(par); end

if ismember('k',flags); PM_makespec_native(par); end

if ismember('z',flags); PM_ArtScansToOns(par); end

if ismember('f',flags); PM_meanFuncs(par); end

if ismember('j',flags); PM_subUtil(par); end

% Modeling, etc
%if ismember('o',flags); par_makemoveonsets(par, (1:3)); end  %EDITED FOR TOP 3 CONF BINS!!
    if ismember('v',flags); PM_MakeRegsLocMVPA(par); end
    if ismember('r',flags); EF_MakeRegs_ON2(par); end
    
%for t = 1:length(par.Tasks)
    if ismember('p',flags); EF_mod_spec(par); end
    if ismember('e',flags); PM_mod_est(par); end
    if ismember('t',flags); EF_setcontrasts(par); end
%end

cd(origdir);


% mail out...
if ismember('M',flags)
     try
        pss = dir([par.logdir filesep '*.ps']);
        [psdts, psord] = sortrows(vertcat(pss.datenum));
        psfile = [par.logdir filesep pss(psord(end)).name];
        setpref('Internet','SMTP_Server',par.smtpserv);
        setpref('Internet','E_mail',par.fromaddy);
        sendmail(par.toaddy, [par.substr ' complete'],...
                ['Complete as of ' datestr(now, 0)], {psfile, parfile});
     catch
        fprintf('Mailing failed');
     end
end
fprintf('\n%s done.\n', par.substr);
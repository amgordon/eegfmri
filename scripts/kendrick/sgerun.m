function sgerun(command,name,wantworkspace,jobindices,priority,flags,memusage)

% function sgerun(-1)
%
% this is equivalent to:
%   sgerun([],[],[],-1);
% which means to use runmatlabinbg.pl (with multiple computational threads) 
% instead of the SGE system.  see below for more details.
%
%   OR
%
% function sgerun(command,name,wantworkspace,jobindices,priority,flags,memusage)
%
% <command> (optional) is a string with the MATLAB code to run.
%   default is [] which means to use an interactive prompt to get the code.
% <name> (optional) is a string with the name of the job, like "job_<name>".
%   default is [] which means to use an interactive prompt to get the name.
%   special case is 0 which means generate a random string.
%   valid names are like ^\w+$.
% <wantworkspace> (optional) is whether to make the current workspace
%   available to the code.  default is [] which means to use an interactive
%   prompt to get the answer.  be careful about saving huge workspaces,
%   as this may be inefficient.
% <jobindices> (optional) has four cases:
%   (1) indices of jobs to do.  indices should be positive integers.  
%       if supplied, we farm out multiple jobs.  the command 
%       that is run by each job is
%         jobindex = X;
%         <command>
%       where X ranges through the elements of <jobindices>.
%       each job gets a name like "job_<name>_X".
%   (2) {A B C D} where
%         A is the name of a variable to use for indexing, like 'ix'
%         B is a string with a vector of indices, like '1:100'
%         C is the number of jobs to create
%         D is a cell vector of variable names, like {'result' 'result2'}
%       in this case, we automatically split up the for-loop implied by
%       A and B.  we assume that the code within the for-loop create the
%       variables named in D.  we also assume that the variables are suitable
%       for merging, as described in loadmulti.m.  (basically, regular matrices
%       should have NaNs for entries that haven't been computed, and cell matrices
%       should have [] for entries that haven't been computed.)  after all jobs 
%       finish, we merge the results from different jobs and then assign variables 
%       in the caller's workspace.  the actual command that is run by each job is 
%       a bit tricky but the user does not need to worry about it.
%   (3) [] which means to do nothing special.
%   (4) -1 which means to use runmatlabinbg.pl (with multiple computational threads) 
%       instead of the SGE system.  in this case, the <priority>, <flags>, and 
%       <memusage> inputs are ignored, and an output file is written to ~/sgeoutput/job_X.out,
%       which is renamed to .out.done upon completion.  make sure to add the path to 
%       runmatlabinbg.pl to your shell file (e.g. .cshrc).
%   default: [].
% <priority> (optional) is an integer in -1024 (lowest) to 1023 (highest).
%   this controls the priority of the job relative to other jobs that you
%   (the user) may have pending.  default: 0.
% <flags> (optional) is a string with additional flags to pass to qsub.
%   this is useful for setting resource requirements for your job.
%   for example, <flags> could be '-l hostname=azure'.  for a list of
%   possible resources, try "qconf -sc".  default: ''.
% <memusage> (optional) is a positive integer indicating the number of megabytes of 
%   memory that you estimate that your job will use.  default: 100.
%
% the purpose of this function is to make it easy to deploy MATLAB code on the
% Sun Grid Engine (SGE).  see also sgestat.m.
%
% basically, what happens is:
%   first, an .m file is created in ~/sgeoutput/ with the code you wish to run, with
%     some various fprintf statements before and after your code.  the file is named 
%     like "job_<name>.m".  
%   then, we use the qsub utility to create an SGE job named like "job_<name>".
%     the actual command associated with the SGE job is "~/matlabsge.sh job_<name>",
%     which, assuming the ~/matlabsge.sh script is setup correctly, deletes any
%     existing output and error files (see below), and then simply runs 
%     the MATLAB command "job_<name>".  outputs and errors from the job are
%     written to the ~/sgeoutput/ directory in files named like "job_<name>.oN"
%     and "job_<name>.eN".  the qsub call that we use is written to a text
%     file "job_<name>.command" --- to run this file, just type something like
%     "sh job_<name>.command".
%   if <wantworkspace>, we save the workspace to the ~/sgeoutput/ directory
%     in a file named like "job_<name>.mat", and ensure that the first thing that
%     happens in the SGE job is the loading of this .mat file.
%
% one of the flags that we use for qsub is a flag that sends e-mail whenever a
% job is rescheduled or aborted (e.g. if an error occurs).  in order to use
% this feature, you must run
%   setpref('kendrick','email','user@host');
% with the e-mail address that you want to use.  if the e-mail address preference is
% not available, we issue a warning and then do not use the feature.
%
% another flag that we use for qsub is -q <queue> where <queue> is the queue to 
% submit the job to.  the default is 'all.q'.  if you want to submit to a different
% queue, you must run
%   setpref('kendrick','sgequeue','blah.q');
% where blah.q is the name of the queue that you want.  
%
% in order to use sgerun.m, you must prepare the following:
%   1. make the directory ~/sgeoutput/
%   2. make the script ~/matlabsge.sh (a sample version is included below)
%   3. permanently add ~/sgeoutput/ to your MATLAB path
%   4. ensure that the SGE utility qsub is available on the current machine.
%      this should involve adding something like this to your .cshrc file:
%        if (-f /usr/share/gridengine/default/common/settings.csh) then
%          source /usr/share/gridengine/default/common/settings.csh
%        endif 
%
% note that you should clean out ~/sgeoutput/ periodically by deleting the various
% files that are associated with jobs that are no longer needed.  otherwise, lots 
% of files will pile up in the directory.
%
% the case when <jobindices> is {A B C D} is a little tricky, so we include here
% an example of that sort of usage:
%
%   here is a simple for-loop done in the normal way.
%
%   result = [];
%   for vx=1:100
%     result(vx) = vx.^2;
%   end
%   isequal(result,(1:100).^2)
%
%   here is the same loop, but dispatched to the SGE.
%
%   result = NaN*zeros(1,100);
%   sgerun([],0,1,{'vx' '1:100' 5 {'result'}});
%     result(vx) = vx.^2;
%   .
%   isequal(result,(1:100).^2)
%
%   in the above call, we first initialize the result using NaNs to indicate entries
%   that have not been computed.  we split the 100 indices to be processed into 5 different jobs.
%   each job works on a different group of the indices.  the 'result' variable from the
%   various jobs are then merged back together and assigned in the current workspace.
%
% sample ~/matlabsge.sh script:
% =======================================
% #!/bin/sh
% 
% # call like matlabsge.sh <command>
% #
% # this will delete any existing output and error files and then run <command> using MATLAB.
% # edit this script to use the particular version of MATLAB that you desire.
% 
% rm -f ~/sgeoutput/$JOB_NAME.o*
% rm -f ~/sgeoutput/$JOB_NAME.e*
% echo "$1" | /usr/local/matlab/r2009b/bin/matlab -nosplash -nodesktop -nodisplay -nojvm -singleCompThread
% =======================================
%
% history:
% 2012/10/24 - to improve speed, no longer make (chmod) the .command files executable;
%              to improve speed, the job submissions issued by this function do not first
%              delete the .o and .e files.  this is unnecessary anyway since this function
%              tries to not allow a job to be submitted for which files in the sgeoutput
%              directory already exist.
% 2011/09/22 - give error if a job name is too long
% 2010/11/19 - add the feature of submitting to a specific queue; implement sgerun(-1) case; no longer writing to recentjobs.txt file
% 2010/09/28 - add -S /bin/sh to the qsub calls

% NOTES:
% maybe we should just have a single chmod call at the end!

% internal constants
masterdir = '~/sgeoutput/';
masterscript = '~/matlabsge.sh';
numletters = 3;      % number of letters in randomly generated job name
maxsize = 10000000;  % warn if workspace is bigger than this number of bytes

% deal with command input
if exist('command','var') && isequal(command,-1)
  command = [];
  name = [];
  wantworkspace = [];
  jobindices = -1;
end
if ~exist('command','var') || isempty(command)
  fprintf('What command do you wish to run?\n');
  command = inputmulti;
end

% deal with name input
if ~exist('name','var') || isempty(name)
  name = [];
end
if isempty(name)
  while 1
    name = input('What is the name of this job? [default is to randomly generate a name]\n','s');
    if isempty(name)
      name = 0;
      break;
    end
    if isempty(regexp(name,'^\w+$'))
      fprintf('Invalid name. Try again.\n');
    elseif exist([masterdir 'job_' name '.m'],'file') || exist([masterdir 'job_' name '_1.m'],'file')
      fprintf('This job name already exists. Try again.\n');
    else
      break;
    end
  end
end
if isequal(name,0)
  isbad = 1;
  while isbad
    name = randomword(numletters);
    isbad = exist([masterdir 'job_' name '.m'],'file') || exist([masterdir 'job_' name '_1.m'],'file');
  end
end
if exist([masterdir 'job_' name '.m'],'file') || exist([masterdir 'job_' name '_1.m'],'file')
  error('The job name already exists.');
end
name = ['job_' name];

% deal with wantworkspace input
if ~exist('wantworkspace','var') || isempty(wantworkspace)
  wantworkspace = input('Do you want the workspace to be available to your job? [0=no [default], 1=yes]\n');
  if isempty(wantworkspace)
    wantworkspace = 0;
  end
  assert(isequal(wantworkspace,0) | isequal(wantworkspace,1));
end

% deal with other inputs
if ~exist('jobindices','var') || isempty(jobindices)
  jobindices = [];
end
if ~exist('priority','var') || isempty(priority)
  priority = 0;
end
if ~exist('flags','var') || isempty(flags)
  flags = '';
end
if ~exist('memusage','var') || isempty(memusage)
  memusage = 100;
end
assert(isint(priority) && priority >= -1024 && priority <= 1023);
assert(isint(memusage) && memusage >= 1);

% % ok. make a file indicating that a job is being created.
% unix_wrapper(sprintf('touch %s',[masterdir name]),0);

% deal with the workspace
if wantworkspace

  % check size of workspace
  a = evalin('caller','whos');
  workspacesize = sum(cat(1,a.bytes));
  
  % give warning if the workspace is big
  if workspacesize > maxsize
    warning('We are saving to disk a workspace that is larger than 10 MB!');
  end

  % save caller's workspace to the special .mat file
  fprintf('saving workspace to disk...');
  evalin('caller',sprintf('save(''%s'');',[masterdir name '.mat']));
  fprintf('done.\n');

  % prepend a load command
  prefix = sprintf(['load(''' masterdir name '.mat''); ']);

else
  prefix = [];
end

% more prefix and suffix stuff
prefix = [prefix 'fprintf(''The job is %s (%s). The current host is %s. The PID is %d. ' ...
                 'The start time is %s.\n'',getenv(''JOB_NAME''),getenv(''JOB_ID''),gethostname,getpid,datestr(now)); '];
suffix = 'fprintf(''The job is %s (%s). The end time is %s.\n'',getenv(''JOB_NAME''),getenv(''JOB_ID''),datestr(now));';

% hack the single-job case to be like the multiple-job case
if isempty(jobindices)
  issinglecase = 1;
  jobindices = 1:1;
elseif isequal(jobindices,-1)
  issinglecase = 2;
  jobindices = 1:1;
else
  issinglecase = 0;
end

% deal with email flag
email = getpref('kendrick','email','');
if isempty(email)
  warning('no e-mail address preference found. please see "help sgerun".');
  emailflag = '';
else
  emailflag = ['-m a -M ' email];
end

% deal with queue name
queuename = getpref('kendrick','sgequeue','all.q');

% finally, let's do it
if iscell(jobindices)

  % give names to variables
  indexvar = jobindices{1};
  indices = jobindices{2};
  numchunks = jobindices{3};
  varstosave = jobindices{4};
  
  % construct the command
  command0 = {prefix sprintf('for %s = chunking(%s,%d,str2double(getenv(''SGE_TASK_ID'')))',indexvar,indices,ceil(length(eval(indices))/numchunks)) command 'end' ...
                             ['save([''' masterdir name '_result''' ' getenv(''SGE_TASK_ID'') ' '''.mat'']' cell2str2(varstosave) ');'] suffix};
  name0 = name;
  if length(name0) > namelengthmax
    error(sprintf('The job name (%s) is too long.',name0));
  end
  
  % write the command to an .m file
  savetext([masterdir name0 '.m'],command0);
  
  % prep
  cmd0 = sprintf('rm -f ~/sgeoutput/%s_*.[oe]*; ',name0);
  cmd1 = sprintf('qsub -t 1-%d -sync y -N %s -o %s -e %s -p %d -l virtual_free=%dM %s %s -S /bin/sh -q %s %s %s; ', ...
                 name0,numchunks,name0,masterdir,masterdir,priority,memusage,emailflag,flags,queuename,masterscript,name0);
  cmd2 = '';%%%%%sprintf('echo "%s" >> ~/sgeoutput/recentjobs.txt; ',cmd1);
  cmd3 = sprintf('echo "%s%s" > ~/sgeoutput/%s.command; ',cmd0,cmd1,name0);     %%%%chmod +x ~/sgeoutput/%s.command;    ,name0
  allcmd = [cmd2 cmd3 cmd1];  % do the echos before the actual command
  
  % do it!
  unix_wrapper(allcmd,1);
  
  % suck in results and assign in the caller's workspace
  files = sortnumerical(matchfiles([masterdir name '_result*.mat']));
  assert(length(files)==numchunks);
  for pp=1:length(varstosave)
    assignin('caller',varstosave{pp},loadmulti(files,varstosave{pp}));
  end

else

  % do it
  allcmd = '';
  cnt = 0;
  for p=jobindices
  
    % prep
    switch issinglecase
    case {1 2}
      command0 = {prefix command suffix};
      name0 = name;
    case 0
      command0 = {prefix sprintf('jobindex = %d;',p) command suffix};
      name0 = sprintf([name '_%d'],p);
    end
    if length(name0) > namelengthmax
      error(sprintf('The job name (%s) is too long.',name0));
    end
  
    % write the command to an .m file
    savetext([masterdir name0 '.m'],command0);
    
    % prep
    cnt = cnt + 1;
    if issinglecase==2
      cmd0 = sprintf('rm -f ~/sgeoutput/%s.[oe]*; ',name0);
      cmd1 = sprintf('/usr/bin/nohup runmatlabinbg.pl 1 %s ~/sgeoutput/%s.out > /dev/null &',name0,name0);
      cmd2 = ''; %%%%%%%sprintf('echo "%s" >> ~/sgeoutput/recentjobs.txt; ',cmd1);
      cmd3 = sprintf('echo "%s" > ~/sgeoutput/%s.command; ',cmd0,cmd1,name0);   %%%%%chmod +x ~/sgeoutput/%s.command;   ,name0
      allcmd = [cmd1 cmd2 cmd3];
      
      % do it?
      unix_wrapper(allcmd,1);
    else
      cmd0 = sprintf('rm -f ~/sgeoutput/%s.[oe]*; ',name0);
      cmd1 = sprintf('qsub -N %s -o %s -e %s -p %d -l virtual_free=%dM %s %s -S /bin/sh -q %s %s %s; ',name0,masterdir,masterdir,priority,memusage,emailflag,flags,queuename,masterscript,name0);
      cmd2 = ''; %%%%%%%sprintf('echo "%s" >> ~/sgeoutput/recentjobs.txt; ',cmd1);
      cmd3 = sprintf('echo "%s%s" > ~/sgeoutput/%s.command; ',cmd0,cmd1,name0);  %%%%%chmod +x ~/sgeoutput/%s.command;    ,name0
      allcmd = [allcmd cmd1 cmd2 cmd3];
      
      % do it?
      if cnt==20 || p==jobindices(end)
        unix_wrapper(allcmd,1);
        allcmd = '';
        cnt = 0;
      end
    end
    
  end
  
  % do it!
  fprintf('\n\nall jobs successfully created!\n');

end

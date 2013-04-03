function runme(varargin)

% function runme
% function runme(mode)
% function runme(mode,command,opt,name)
% function runme(mode,indices,command,opt,name)
%
% <mode> (optional) is
%   0 means go ahead and issue the "runmatlabinbg.pl" command, like:
%     /usr/bin/nohup runmatlabinbg.pl 1 job_userid_<name> /home/knk/jobs/job_userid_<name>.out > /dev/null &
%   1 means let the queue daemon deal with running the job, by creating a file:
%     /home/knk/queue/job_userid_<name>
%   default: 0.  note that userid is determined by running whoami.
% <command> (optional) is a string with the command to run.
%   default is [] which means to use interactive prompt.
% <opt> (optional) is some data that the command requires.
%   if supplied, we save this data as 'opt' and make sure 'opt'
%   is available to the command when run.
%   default is [] which means to do nothing special.
% <name> (optional) is a string.
%   default is [] which means to use interactive prompt.
%   special case is 0 which means generate a random string.
% <indices> (optional) is the indices of jobs to do.  if supplied, 
%   we farm out multiple jobs.  the command that is run by each job is
%   determined by sprintf(command,x) where x ranges through <indices>,
%   and each job gets a name like "<name>X" where X is an element of <indices>.
%   default is [], which means do nothing special.  
%
% we need to document more about the queueing setup!!!
%
% make sure to add the path to runmatlabinbg.pl to your shell file (e.g. .cshrc).
%
% see also queuedaemon.m.
%
% history:
% 2010/08/11 - (big change; see also queuedaemon.m) input format totally changed.  so calls have to be modified accordingly.

% NOTES FOR FUTURE DOCUMENTATION:
% /home/knk/jobs/%s.m
% /home/knk/jobs/%s.mat

% input
switch nargin
case 0
  mode = [];
  command = [];
  opt = [];
  name = [];
  indices = [];
case 1
  mode = varargin{1};
  command = [];
  opt = [];
  name = [];
  indices = [];
case 4
  mode = varargin{1};
  command = varargin{2};
  opt = varargin{3};
  name = varargin{4};
  indices = [];
case 5
  mode = varargin{1};
  indices = varargin{2};
  command = varargin{3};
  opt = varargin{4};
  name = varargin{5};
end
if isempty(mode)
  mode = 0;
end

% read input
if isempty(command)
  fprintf('what command do you wish to run?\n');
  command = inputmulti;
end

% get filename
if isempty(name)
  while 1
    id0 = input('what is the name of this job?\n','s');
    if isempty(regexp(id0,'^\w+$'))
      fprintf('invalid name. try again.\n');
    else
      break;
    end
  end
elseif isequal(name,0)
  id0 = randomword(10);
else
  id0 = name;
end

% figure out user id
userid = matchword(unix_wrapper('whoami',0));

% single-job case
if isempty(indices)
  runme_helper(mode,command,sprintf('job_%s_%s',userid,id0),opt);

% multiple-jobs case
else
  for x=indices
    runme_helper(mode,sprintf(command,x),sprintf('job_%s_%s_%d',userid,id0,x),opt);
  end
end

%%%%%

function runme_helper(mode,command,id,opt)

% construct filenames and check for existence
file = sprintf('/home/knk/jobs/%s.m',id);
fileout = sprintf('/home/knk/jobs/%s.out',id);
fileoutdone = sprintf('/home/knk/jobs/%s.out.done',id);
filequeue = sprintf('/home/knk/queue/%s',id);
filemat = sprintf('/home/knk/jobs/%s.mat',id);
if exist(file,'file') || exist(fileout,'file') || exist(fileoutdone,'file') || exist(filequeue,'file') || exist(filemat,'file')
  error('name already exists');
end

% deal with opt
if ~isempty(opt)
  save(filemat,'opt');
  command = ['load(''' filemat '''); ' command];
end

% write input to file
savetext(file,command);

% do it
switch mode
case 0

  % run matlab in bg
  [status,result] = unix(sprintf('/usr/bin/nohup runmatlabinbg.pl 1 %s %s > /dev/null &',id,fileout));
  assert(status==0);

case 1
  
  % make a dummy placeholder file for the queue daemon
  savetext(filequeue,id);

end

% report
fprintf('%s created.\n',id);

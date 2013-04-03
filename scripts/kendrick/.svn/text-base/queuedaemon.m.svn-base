function queuedaemon

% function queuedaemon
%
% run indefinitely.  use /home/knk/queuestate.mat if it exists; otherwise, use initializations.
% stop if a file named /home/knk/queuestop exists.  record state in /home/knk/queuestate.mat.
%
% commands to run are indicated by the existence of files named "/home/knk/queue/job_userid_xyz".
% after we detect the existence of files, we dispatch the oldest command
%   and then delete the corresponding file "/home/knk/queue/job_userid_xyz".
%   to dispatch commands, we ssh to a machine logging in as userid and use
%   runmatlabinbg.pl to run the command.  see code for the specific call that we use.
% the output of a given command is written to "/home/knk/jobs/job_userid_xyz.out".
% when a command is finished, the output is renamed to "/home/knk/jobs/job_userid_xyz.out.done".
%
% run like this:
%   /usr/bin/nohup runmatlabinbg.pl 0 queuedaemon /dev/null > /dev/null &
%
% in order to get the ssh mechanism working, you have to make it so that passwords are not prompted:
%   on local machine:
%     cd .ssh
%     ssh-keygen -t dsa
%     (leave password blank)
%     scp id_dsa.pub user@remote:~/.ssh/
%   on remote machine:
%     cd .ssh
%     cat id_dsa.pub >> authorized_keys2
%     chmod 640 authorized_keys2
%     rm id_dsa.pub
%
% also, to run jobs, you have to add /home/knk/jobs to your matlab path.
%
% see also runme.m.
%
% history:
% 2010/08/11 - (big change; see also runme.m) require job_userid_xyz format.  when dispatching, ssh using userid.  remove local running option.

% internal constants
statefile = '/home/knk/queuestate.mat';
stopfile = '/home/knk/queuestop';

% if state file doesn't exist, use some inits
if ~exist(statefile,'file')
  % constants
  maxproc = [4 32 0 0];  % maximum queue jobs to run at a time  % tan,purple support 8
  machine = {'chroma' 'azure' 'tan' 'purple'};
  
  % init
  cnt = [0 0 0 0];   % number of queue jobs running right now
  running = {};  % record of the queue commands that we're running
  runningmachine = [];  % record of which machine it's running on

% use state if it exists
else
  load(statefile);
end

% prep
warning off;

% do it
while 1

  % save state
  save(statefile);
  
  % stop queue?
  if exist(stopfile,'file')
    delete(stopfile);
    return;
  end

  % wait a little
  pause(5);

  % see if any jobs have completed, and adjust cnt and running accordingly
  isdone = zeros(1,length(running));
  for p=1:length(running)
    id = running{p};
    fileoutdone = sprintf('/home/knk/jobs/%s.out.done',id);
    isdone(p) = exist(fileoutdone,'file');
    if isdone(p)
      %fprintf('%s is DONE.\n',id);
    end
  end
  ix = find(isdone~=0);  % if not 0, you exist, and therefore done
  for p=1:length(ix)
    cnt(runningmachine(ix(p))) = cnt(runningmachine(ix(p))) - 1;
  end
  running(ix) = [];
  runningmachine(ix) = [];

  if ~isempty(ix)
    %fprintf('after adjusting for done, cnt=%s.\n',mat2str(cnt));
    %fprintf('after adjusting for done, runningmachine=%s.\n',mat2str(runningmachine));
  end
  
  % if we can add jobs, let's do it
  if any(cnt < maxproc)
    %fprintf('1');

    % see if there are any jobs to run
    files = matchfiles('/home/knk/queue/*','tr');
    if ~isempty(files)
      
      % get the first one
      filequeue = files{1};
      frags = strsplit(filequeue,'/');
      id = frags{end};
      userid = firstelc(firstelc(regexp(id,'job_(.+?)_','tokens')));
      fileout = sprintf('/home/knk/jobs/%s.out',id);
      fileoutdone = sprintf('/home/knk/jobs/%s.out.done',id);
      
      for p=randperm(length(maxproc))  % check in a random order!
        if cnt(p) < maxproc(p)
          % check load
          [status,result] = unix(sprintf('ssh -x -o ConnectionAttempts=10 %s "cat /proc/loadavg"',machine{p}));
          if status~=0
            continue;
          end
          temp = str2num(result);
          if temp(1) > maxproc(p)
            continue;
          end

          % run matlab in bg remotely
            % -o StrictHostKeyChecking=no -o PasswordAuthentication=no 
          %fprintf('RUN %s remotely...',id);
          [status,result] = unix(sprintf('ssh -x -o ConnectionAttempts=10 -f %s@%s "runmatlabinbg.pl 0 %s %s"',userid,machine{p},id,fileout));
          if status~=0
            continue;  % hm, didn't work, so let's act like we didn't try it
          end
          %fprintf('ok.\n');
          
          % if .out.done already exists for some reason, delete it
          delete(fileoutdone);
          
          % update
          cnt(p) = cnt(p) + 1;
          running{end+1} = id;
          runningmachine(end+1) = p;
          %fprintf('after running, cnt=%s.\n',mat2str(cnt));
          %fprintf('after running, runningmachine=%s.\n',mat2str(runningmachine));

          % remove from queue directory
          delete(filequeue);
          
          % we're done
          break;
        end
      end

    end

  end

end




%           if isequal(machine{p},'')  % THIS IS OBSOLETE SINCE WE JUST SSH TO THE HOST MACHINE
%             % check load
%             [status,result] = unix('cat /proc/loadavg');
%             assert(status==0);
%             temp = str2num(result);
%             if temp(1) > maxproc(p)
%               continue;
%             end
%           
%             % run matlab in bg locally
%             %fprintf('RUN %s locally...',id);
%             [status,result] = unix(sprintf('/usr/bin/nohup runmatlabinbg.pl 0 %s %s > /dev/null &',id,fileout));
%             assert(status==0);
%             %fprintf('ok.\n');
%           else

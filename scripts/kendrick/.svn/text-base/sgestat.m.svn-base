function sgestat(flag)

% function sgestat(flag)
%
% <flag> (optional) is a non-negative integer with bit
%   1 indicating whether to show finished jobs
%   2 indicating whether to show running jobs
%   3 indicating whether to show unstarted jobs
%   default: 7.
%
% print out information pertaining to all known SGE jobs.
% see sgerun.m for more details.
%
% the printed information should be self-explanatory, except for "last update":
% in this case, there are two durations.  the first is the duration between
% when the job was started and the last modification of the .o file.
% the second is the duration between the last modification and the 
% current time, which is useful for detecting jobs that are not updating
% for some weird reason.

% internal constants
formatline = '% 50.50s  % 6.6s  % 23.23s  % 23.23s  % 6.6s  % 15.15s  % 6.6s  % 32.32s\n';
headernum = 50;

% input
if ~exist('flag','var') || isempty(flag)
  flag = 7;
end

% all the .m, .o, and .e files that exist
allfilesM = dir('~/sgeoutput/*.m'); 
allfilesO = dir('~/sgeoutput/*.o*'); 
allfilesE = dir('~/sgeoutput/*.e*'); 
  %        name: 'job_adc_1.e1696'
  %        date: '17-Sep-2010 23:43:34'    #!
  %       bytes: 0
  %       isdir: 0
  %     datenum: 734398.988587963
assert(length(allfilesO)==length(allfilesE),'number of .o files does not equal number of .e files');

% all the job names
jobnames = arrayfun(@(x) firstelc(regexp(x.name,'(job.+)\.m','tokens')),allfilesM);  % {'job_rk2dc_9' ...}

% all the jobs that have been started and their job numbers and which ones have non-empty error files
startednames = arrayfun(@(x) firstelc(regexp(x.name,'(job.+)\.','tokens')),allfilesO);  % {'job_rk2dc_9' ...}
startedjobnums = arrayfun(@(x) firstelc(regexp(x.name,'job.+\.o(.+)$','tokens')),allfilesO);  % {'1965' ...}
startedhaserror = cat(2,allfilesE.bytes)~=0;  % [0 0 0 0 1 0 0 ...]
startedmods = {allfilesO.date};  % {'17-Sep-2010 23:43:34' ...}
if isempty(startednames)
  startednames = {};
end

% sort
[jobnames,ix] = sortnumerical(jobnames);
[startednames,ix] = sortnumerical(startednames); startedjobnums = startedjobnums(ix); startedhaserror = startedhaserror(ix); startedmods = startedmods(ix);

% load in information from the .o files
if length(dir('~/sgeoutput/*.o*'))==0
  headinfo = {};
  tailinfo = {};
else
  headinfo = regexp(unix_wrapper('head -100 ~/sgeoutput/*.o* | grep "The current host is"',0,0), ...
                     '.*?The job is (.+?) \((.+?)\)\. The current host is (.+?)\. The PID is (\d+?)\. The start time is (.+?)\.','tokens');
    % {{'job_flsdfkj_1' '988' 'slate.stanford.edu' '27934' '18-Sep-2010 05:32:56'} ...}
  tailinfo = regexp(unix_wrapper('tail -n 10 ~/sgeoutput/*.o* | grep "The end time is"',0,0), ...
                   '.*?The job is (.+?) \((.+?)\)\. The end time is (.+?)\.','tokens');
    % {{'job_flsdfkj_1' '988' '18-Sep-2010 05:32:56'} ...}
end
headnames = cellfun(@(x) x(1),headinfo);
tailnames = cellfun(@(x) x(1),tailinfo);
if isempty(headnames)
  headnames = {};
end
if isempty(tailnames)
  tailnames = {};
end
[headnames,ix] = sortnumerical(headnames); headinfo = headinfo(ix);
[tailnames,ix] = sortnumerical(tailnames); tailinfo = tailinfo(ix);

% save header line
jobname = 'job name';
sgeid = 'SGE ID';
whenstarted = 'when started';
host = 'host';
pid = 'PID';
lastupdate = 'last update';
haserror = 'error?';
whenfinished = 'when finished';
header = sprintf(formatline,jobname,sgeid,whenstarted,host,pid,lastupdate,haserror,whenfinished);

% DO IT
cnt = 1;
nfinished = 0;
nrunning = 0;
nunstarted = 0;
for p=1:length(jobnames)

  % calc
  startix = find(ismember(startednames,{jobnames{p}}));
  headix = find(ismember(headnames,{jobnames{p}}));
  tailix = find(ismember(tailnames,{jobnames{p}}));

  % figure out values
  jobname = jobnames{p};
  if isempty(headix)
    whenstarted = '';
    host = '';
    pid = '';
  else
    whenstarted = headinfo{headix}{5};
    host = headinfo{headix}{3};
    pid = headinfo{headix}{4};
  end
  if isempty(startix)
    sgeid = '';
    lastupdate = '';
    haserror = '';
  else
    sgeid = startedjobnums{startix};
    lastupdate = startedmods{startix};
    haserror = choose(startedhaserror(startix),'ERROR','');

    if isempty(headix)
      nhours = 0;
      nmins = 0;
    else
      hh = (datenum(lastupdate)-datenum(whenstarted))*24;
      nhours = floor(hh);
      nmins = round((hh-floor(hh))*60);
    end

    hh = (datenum(now)-datenum(lastupdate))*24*60;
    nminsB = round(hh);

    lastupdate = sprintf('%dh%dm + %dm',nhours,nmins,nminsB);
    
    if ~isempty(tailix)
      lastupdate = '-----';
    end

  end
  if isempty(tailix)
    whenfinished = '';
  else
    whenfinished = tailinfo{tailix}{3};

    hh = (datenum(whenfinished)-datenum(whenstarted))*24;
    nhours = floor(hh);
    nmins = round((hh-floor(hh))*60);

    whenfinished = [whenfinished sprintf(' (%dh%dm)',nhours,nmins)];
  end
  
  % increment counters
  nfinished = nfinished + ~isempty(tailix);
  nrunning = nrunning + (isempty(tailix) && ~isempty(startix));
  nunstarted = nunstarted + isempty(startix);
  
  % if finished job and the first bit is 0
  if ~isempty(tailix) && bitget(flag,1)==0
    continue;
  end
  % if running job and the third bit is 0
  if isempty(tailix) && ~isempty(startix) && bitget(flag,2)==0
    continue;
  end
  % if unstarted job and the second bit is 0
  if isempty(startix) && bitget(flag,3)==0
    continue;
  end
  
  % print out a header every N lines
  if mod(cnt,headernum)==1
    fprintf([repmat('=',1,length(header)) '\n']);
    fprintf(header);
    fprintf([repmat('=',1,length(header)) '\n']);
  end

  % print it
  fprintf(formatline,jobname,sgeid,whenstarted,host,pid,lastupdate,haserror,whenfinished);
  cnt = cnt + 1;
  
end

% final stats
fprintf('There are %d finished jobs.\n',nfinished);
fprintf('There are %d running jobs.\n',nrunning);
fprintf('There are %d unstarted jobs.\n',nunstarted);

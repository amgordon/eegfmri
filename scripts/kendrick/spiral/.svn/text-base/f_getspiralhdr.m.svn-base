function [nslc,necho,nviews,nfrms,npts,rhptsize,rec_start,rec_end,rhuser] = ...
        f_getspiralhdr(raw_file_name)
% M-file function to get header information from Atsushi's spiral PSD
% [nslc,necho,nviews,nfrms,npts,rhptsize,rec_start,rec_end,rhuser]=f_getspiralhdr(raw_file_name)
% rhuser0 is in rhuser(1) and so on (shifted by one)

fid=fopen(raw_file_name,'r','ieee-le');
if fid == -1, error('File read Error'), end
% Get file size information from header
revision = fread(fid,1,'float');
run_int = fread(fid,1,'int');
scan_int = fread(fid,1,'short');
run_ascii = sprintf('%s',fread(fid,6,'char'));
date_ascii = sprintf('%s',fread(fid,10,'char'));
time_ascii = sprintf('%s',fread(fid,8,'char'));
logo_ascii = sprintf('%s',fread(fid,10,'char'));

clear hdr
status=fseek(fid,0,'bof');
hdr = fread(fid, 39984/2, 'short');
rhptsize = hdr(42);
%%% if (rhptsize > 4) % We must be byte swapped... reopen as ieee-le 
if (sscanf(run_ascii,'%d')==run_int) % We must be byte swapped... reopen as ieee-le 
   fclose(fid);
   fid=fopen(raw_file_name,'r','ieee-le');
   disp('Reopening file as ieee-le')
   revision = fread(fid,1,'float');
   run_int = fread(fid,1,'int');
   scan_int = fread(fid,1,'short');
   run_ascii = fread(fid,6,'char')';
   if fid == -1, error('File read Error'), end
% Get file size information from header
   clear hdr
   status=fseek(fid,0,'bof');
   hdr = fread(fid, 39984/2, 'short');
   rhptsize = hdr(42);
   if (rhptsize > 4)
      error('cannot figure out this byte swapping stuff');
   end
end

disp(sprintf('Header Revision %g',revision))
disp(sprintf('Run Number (Long) %d, (ASCII) %s, Pfilename %s', ...
              run_int,         run_ascii, raw_file_name))
disp(sprintf('Scan Number (Integer) %d',scan_int))
disp(sprintf('Date & Time : %s %s',date_ascii,time_ascii))
disp(sprintf('Logo        : %s',logo_ascii))
no_proc = rem(hdr(25),2);, if(no_proc == 1), fprintf('NO_PROC data\n');, end;
nslc = hdr(35);
necho = hdr(36);
nviews = hdr(38);
npts = hdr(41);

%  fprintf('npts = %d   nviews = %d   ',npts, nviews);
%  fprintf('nslc = %d   necho = %d   ',nslc,necho);
%  fprintf('npts = %d   nviews = %d   ',npts, nviews);
%  fprintf('nslc = %d   necho = %d  \n',nslc,necho);
disp('Find receiver numbers (Informational; not implemented yet)')
rhdab0s = hdr(101);
rhdab0e = hdr(102);
rhdab1s = hdr(103);
rhdab1e = hdr(104);
rhdab2s = hdr(105);
rhdab2e = hdr(106);
rhdab3s = hdr(107);
rhdab3e = hdr(108);
fprintf('rhdab0s = %d rhdab0e = %d\n',rhdab0s, rhdab0e)
fprintf('rhdab1s = %d rhdab1e = %d\n',rhdab1s, rhdab1e)
fprintf('rhdab2s = %d rhdab2e = %d\n',rhdab2s, rhdab2e)
fprintf('rhdab3s = %d rhdab3e = %d\n',rhdab3s, rhdab3e)

nppv = 2*npts;
nppimg = nppv*(nviews+1);

nppv = 2*npts;
nppimg = nppv*(nviews+1);
% Read the floats
status=fseek(fid,0,'bof');
hdrf = fread(fid, 39984/4, 'float');
rhuser(1:20)=hdrf(55:74);       %rhuser0 to rhuser19
rhuser(21:48)=hdrf(251:278);    %rhuser20 to rhuser48
gdeltat=rhuser(11);             % rhuser10
GAM=rhuser(12);                 % rhuser11
f_offset=rhuser(13);            % rhuser12
tsp=rhuser(14);                 % rhuser13
opmaptype=rhuser(15);           % rhuser14
map_delat=rhuser(16);           % rhuser15
k_max=rhuser(17);               % rhuser16
g_max=rhuser(18)*10;            % rhuser17 in mT/m
s_max=rhuser(19);               % rhuser18
nturns=rhuser(20);              % rhuser19
f  = hdrf (67);                 % frequency offset  rhuser12
dt = hdrf (68)*1e-6;            % sampling time     rhuser13
dx = rhuser(23);                % rhuser22
dy = rhuser(24);                % rhuser23
fprintf('f=%g Hz   dt=%g usec\n',f,dt*1e6);
fprintf('Spiral: k_max = %g, g_max= %g, s_max= %g\n',k_max,g_max,s_max);
fprintf('        deltat=%g turns= %g, GAMMA=%g\n',gdeltat,nturns,GAM)

fprintf('Slice offset = %g, %g\n',dx,dy);
nfrms = rhuser(1);

rec_start=rhdab0s;
rec_end=rhdab0e;

status=fclose(fid);
if status ~= 0, error('File read Error'), end

function greconsevaluateb0(grecons,files,figuredir,n,s,vals,numframes,isinout,gain)

% function greconsevaluateb0(grecons,files,figuredir,n,s,vals,numframes,isinout,gain)
% 
% <grecons> is a string indicating the grecons utility that you wish to use.
%   the following versions have been tested: grecons14 rev 95.
%   no guarantees if you use a different version of grecons.
% <files> is a file pattern (see matchfiles.m) matching one or more .7 files.
%   each .7 file should correspond to one functional run.
% <figuredir> is the directory to write figures to.  this directory should end in
%   an actual directory name because we automatically extract the directory name D
%   and use the prefix 'D_' for figure files.
% <n> is the power of two output matrix size (e.g. 64).  must be greater than 1.
% <s> is the number of slices.  must be greater than 1.
% <vals> is a vector of b0 values to try.  note that the final time-series written 
%   to the .mag file corresponds to the first entry in <vals>.
% <numframes> is
%   [] means normal operation (reconstruct all volumes)
%   N where N is a positive integer means to call grecons with "-f N"
% <isinout> (optional) is whether the data are spiral in-out.  default: 0.
% <gain> (optional) is a number for image scaling (passed to grecons).  default: 1.
%
% reconstruct each run at each b0 value.  for each reconstruction, take the last volume
% returned by grecons (controlled by <numframes>).  write these volumes to run001_001.png,
% etc., where the first number indicates the run number and the second number indicates the
% index of the b0 value.  the same scaling is used for all volumes: [0,X] where X is the 
% 99th percentile of all pixels in all returned volumes.  also, the time-series data for 
% the first specified b0 value (in <vals>) are written to .mag files (to the same 
% directories that contain the original .7 files).  before exiting, we save variables 
% (excluding big ones) to record.mat.  this file is for development purposes only.
%
% please see the code for the exact flags used in the grecons utility.
%
% note that grecons creates files in the directory it is called from as well as in
% the directories that contain the .7 files.  be sure that you do not have two copies 
% of grecons that are trying to reconstruct the same .7 files.  also, we make no attempt
% to clean up the various intermediary files.
%
% some notes on how the case of spiral-in-out data (<isinout>) is handled:
% - for the inspection volumes, we write out two sets of volumes, first the spiral-in 
%   volumes and then the spiral-out volumes.
% - the .mag files will have two sets of volumes, first the spiral-in volumes
%   and then the spiral-out volumes (if <numframes> is specified, there will be 
%   2*<numframes> volumes.)
%
% see also greconsmultifreq.m.
%
% history:
% 2010/06/05 - now save record .mat file
% 2010/05/29 - now we save .mag files corresponding to the first <vals> entry
% 2010/05/26 - initial version.  2:15 pm.

% internal constants
dformat = 'int16';                       % the data format used by grecons
ptile = 99;                              % percentile for scaling purposes
oname = 'run%03d_%03d';                  % file name
rname = 'record';                        % the .mat file we write out

% inputs
if ~exist('isinout','var') || isempty(isinout)
  isinout = 0;
end
if ~exist('gain','var') || isempty(gain)
  gain = 1;
end

% calc
prefix = [stripfile(figuredir,1) '_'];

% make figuredir if necessary
mkdirquiet(figuredir);

% match the filenames
files = matchfiles(files);
assert(~isempty(files),'<files> does not match at least one file');
fprintf('we are processing these files:\n'); fprintf(' %s\n',files{:});

% do it
inspections = zeros(0,0,dformat);  % X x Y x Z x VAL x FILE
final = zeros(0,0,dformat);        % X x Y x time x Z
for z=1:length(files)
  for p=1:length(vals)
    for q=1:s
    
      % do the reconstruction
      unix_wrapper(sprintf('%s %s -s %d -e %d -B 0 -o %.10f -i 1 -G %.10f -N %d -O %s', ...
        grecons,choose(isempty(numframes),'',sprintf('-f %d',numframes)),q,q, ...
        vals(p),gain,n,files{z}));
      a = loadbinary([files{z} '.mag'],dformat,[n n 0]);  % X x Y x time

      % figure out inspections
      if isinout
        if isempty(numframes)
          inspections(:,:,q,[p length(vals)+p],z) = reshape(a(:,:,[end/2 end]),[size(a,1) size(a,2) 1 2]);
        else
          inspections(:,:,q,[p length(vals)+p],z) = reshape(a(:,:,[numframes end]),[size(a,1) size(a,2) 1 2]);
        end
      else
        inspections(:,:,q,p,z) = a(:,:,end);
      end
      
      % figure out final
      if p==1
        if isinout
          if isempty(numframes)
            final(:,:,:,q) = a;
          else
            final(:,:,:,q) = a(:,:,[1:numframes end-numframes+1:end]);  % ignore middle zeros
          end
        else
          final(:,:,:,q) = a;
        end
      end

    end
  end
  
  % save data to the binary file
  savebinary([files{z} '.mag'],dformat,final);

end

% write inspection figures out
mx = prctile(double(inspections(:)),ptile);
for z=1:length(files)
  for p=1:size(inspections,4)  % use size of inspections because it might be double the normal number
    imwrite(uint8(255*makeimagestack(double(inspections(:,:,:,p,z)),[0 mx])),gray(256),sprintf([figuredir '/' prefix oname '.png'],z,p));
  end
end

% save record
clear a final;
save([figuredir '/' prefix rname '.mat']);

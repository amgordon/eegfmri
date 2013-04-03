function adjustgrecons(grecons,files,figuredir,n,s,wantavg,cutoff,mode)

% function adjustgrecons(grecons,files,figuredir,n,s,wantavg,cutoff,mode)
% 
% <grecons> is a string indicating the grecons utility that you wish to use.
%   the following versions have been tested: grecons14 rev 85 and rev 95.
%   no guarantees if you use a different version of grecons.
% <files> is a file pattern matching one or more .7 files (see matchfiles.m).
%   each .7 file should correspond to one functional run.
% <figuredir> is the directory to write figures to.  this directory should end in
%   an actual directory name because we automatically extract the directory name D
%   and use the prefix 'D_' for figure files.
% <n> is the power of two output matrix size (e.g. 64)
% <s> is the number of slices
% <wantavg> is a 3-element vector of 0/1 indicating whether to average 
%   each of the three parameters (foff, wx, wy) across the functional runs.
%   for example, [1 0 0] indicates to average the foff parameter but to 
%   not average the wx and wy parameters.  can be a scalar in which case we apply
%   that value for all three parameters.
% <cutoff> is a 3-element vector indicating the number of cycles per field-of-view 
%   in the slice dimension to use for the low-pass filter cutoff for each of the
%   three parameters (foff, wx, wy).  can be Inf (which in effect means to 
%   perform no smoothing).  can be a scalar in which case we apply that value for 
%   all three parameters.
% <mode> is
%   0 means get parameters, average them (if requested) and smooth them (if requested),
%     and write an inspection figure.  this mode is useful for playing with <cutoff>.
%   1 means also calculate the fieldmaps and brains and write them to figures.
%   2 means also perform the reconstructions and write the 6th and last volumes of
%     each functional run to figures.
%
% in <mode> 0, 1, and 2, we use grecons to determine the three parameters
% foff, wx, and wy, which are defined for each slice and each functional run.
% (we obtain these values by grep and regexp on the verbose output of grecons.)
% we then calculate the desired parameter values by taking the mean across runs
% (if requested) and then smoothing across slices (if requested).  smoothing
% is performed in the space domain (assuming replicated end-points) using a 
% Butterworth filter of order 5 and a certain <cutoff>.  (the filter is constrained 
% to sum to 1 so that DC is not distorted.)  figures summarizing these steps are 
% written to parameter1.png, parameter2.png, and parameter3.png (corresponding to 
% foff, wx, and wy).
%
% in <mode> 1 and 2, we then use grecons to calculate the fieldmaps and
% the brains associated with these fieldmaps.  we write these out for inspection
% purposes to fieldmap%03d.png and brain%03d.png where the index represents
% run number (which corresponds to the alphabetical order of <files>).  the
% same scaling is used for all fieldmaps: [-X,X] where X is the 99th percentile
% of the absolute value of all pixels in all fieldmaps.  the same scaling is used
% for all brains: [0,X] where X is the 99th percentile of all pixels in all brains.
%
% in <mode> 2, we then reconstruct the functional runs, forcing each reconstruction
% to use the fieldmap-correction parameters determined above.  the end result 
% is that we write a .mag file for each functional run (written to the 
% same directory that contains the original .7 file).  we also write out the 6th and
% last volumes for inspection purposes to volume%03da.png and volume%03db.png where
% the index represents the run number.  the same scaling is used for all volumes:
% [0,X] where X is the 99th percentile of all pixels in all of the volumes that we
% are writing out.  (if there are less than six volumes in a functional run, we 
% write out the first and last volumes of that run.)
%
% note that intermediate files will be written to the directory containing
% the .7 files.  we make no guarantees on what these files are.
%
% please see the code for the exact flags used in the grecons utility.
%
% also, it appears that grecons makes some files in the directory from which it is
% called, so beware of some junk files lying around.
%
% see also greconsmultifreq.m.
%
% ALTHOUGH THIS FUNCTION IS "CURRENT" AND WORKS, THIS FUNCTION IS DEPRECATED 
% SINCE IT SERVES NO REAL FUNCTIONALITY.  USE GRECONSMULTIFREQ.M INSTEAD.
%
% NOTE: note that spiral in-out is not handled by this code.
%
% history:
% 2010/05/26 - start tracking.  2:15 pm.

% internal constants
dformat = 'int16';            % the data format used by grecons
fname = 'fieldmap%03d';       % file names for the fieldmaps
bname = 'brain%03d';          % file names for the brains
iname = 'volume%03d';         % file names for the inspections
filterorder = 5;              % order of the Butterworth filter to use
pnames = {'foff' 'wx' 'wy'};  % names of the various parameters
nparams = length(pnames);     % total number of parameters
fvol = 6;                     % the "first" volume to inspect
ptile = 99;                   % percentile for scaling purposes
rname = 'record';             % the .mat file we write out

% inputs
if length(wantavg)==1
  wantavg = repmat(wantavg,[1 nparams]);
end
if length(cutoff)==1
  cutoff = repmat(cutoff,[1 nparams]);
end

% calc
prefix = [stripfile(figuredir,1) '_'];

% make figuredir if necessary
mkdirquiet(figuredir);

% match the filenames
files = matchfiles(files);
assert(~isempty(files),'<files> does not match at least one file');
fprintf('we are processing these files:\n'); fprintf(' %s\n',files{:});

% get fieldmap parameters
params = zeros(s,nparams,length(files));
for p=1:length(files)
  for q=1:s
    result = unix_wrapper(sprintf('%s -s %d -e %d -v -B 2 -i 1 -N %d -O %s | grep fof',grecons,q,q,n,files{p}));
    tokens = regexp(result,'foff =(?<foff>.+?) Hz.+?wx, wy =(?<wx>.+?) (?<wy>.+?) rad','names');
    params(q,:,p) = [str2num(tokens.foff) str2num(tokens.wx) str2num(tokens.wy)];
  end
end

% figure out desired parameters
dparams = zeros(s,nparams,length(files));
for p=1:length(files)
  for q=1:nparams
    if wantavg(q)
      dparams(:,q,p) = tsfilter(mean(params(:,q,:),3)',constructbutterfilter1D(s,cutoff(q),filterorder),[1 0])';
    else
      dparams(:,q,p) = tsfilter(params(:,q,p)',constructbutterfilter1D(s,cutoff(q),filterorder),[1 0])';
    end
  end
end

% inspect parameters
for p=1:nparams
  figureprep; hold on;
  h = plot(squish(params(:,p,:),2));
  h2 = plot(squish(dparams(:,p,:),2),'LineWidth',3);
  legend(h,cellfun(@(x) stripfile(x,1),files,'UniformOutput',0));
  xlabel('slice number'); ylabel('parameter value');
  title(sprintf('parameter %s',pnames{p}));
  figurewrite(sprintf([prefix 'parameter%d'],p),[],[],figuredir);
end

if mode >= 1

  % calculate and load in fieldmaps and brains
  fieldmaps = zeros(n,n,s,length(files));
  brains = zeros(n,n,s,length(files));
  for p=1:length(files)
    unix_wrapper(sprintf('%s -B 2 -i 1 -N %d -O %s',grecons,n,files{p}));
    for q=1:s
      fieldmaps(:,:,q,p) = loadbinary(sprintf([files{p} '.B0.%03d'],q),dformat,[n n]);
      brains(:,:,q,p) = loadbinary(sprintf([files{p} '.%03d'],q),dformat,[n n]);
    end
  end
  
  % write fieldmaps and brains to .png files for inspection
  fmx = prctile(abs(fieldmaps(:)),ptile);
  bmx = prctile(brains(:),ptile);
  for p=1:length(files)
    imwrite(uint8(255*makeimagestack(fieldmaps(:,:,:,p),[-fmx fmx])),sprintf([figuredir '/' prefix fname '.png'],p));
    imwrite(uint8(255*makeimagestack(brains(:,:,:,p),[0 bmx])),sprintf([figuredir '/' prefix bname '.png'],p));
  end
  
  if mode >= 2

    % do the reconstruction (saving inspection volumes)
    inspections = [];
    for p=1:length(files)
    
      % reconstruct each slice, collect them up, save it to a binary file
      vol = zeros(0,0,dformat);
      for q=1:s
        unix_wrapper(sprintf('%s -s %d -e %d -B 1 -o %.10f -i 1 -N %d -O -x %.10f -y %.10f %s',grecons,q,q, ...
          dparams(q,1,p)-params(q,1,p),n,dparams(q,2,p)-params(q,2,p),dparams(q,3,p)-params(q,3,p),files{p}));
        vol(:,:,:,q) = loadbinary([files{p} '.mag'],dformat,[n n 0]);
      end
      savebinary([files{p} '.mag'],dformat,vol);
      
      % collect the inspection volumes
      fvol0 = choose(size(vol,3) < fvol,1,fvol);
      inspections = cat(4,inspections,permute(double(vol(:,:,[fvol0 end],:)),[1 2 4 3]));

    end
      
    % write inspection figures
    bmx = prctile(inspections(:),ptile);
    for p=1:length(files)
      imwrite(uint8(255*makeimagestack(inspections(:,:,:,(p-1)*2+1),[0 bmx])),sprintf([figuredir '/' prefix iname 'a.png'],p));
      imwrite(uint8(255*makeimagestack(inspections(:,:,:,(p-1)*2+2),[0 bmx])),sprintf([figuredir '/' prefix iname 'b.png'],p));
    end

  end

end

% save record
clear vol;
save([figuredir '/' prefix rname '.mat']);





% This is the documentation returned by grecons14 rev 85:
%
% Usage: grecons14 [-s isl1] [-e isl2] [-f nframes] [-B b0map] [-o freq_offset]
%        [-i imgtyp] [-N npix] [-g gamma] [-G gain] [-z zoomf] [-F fermif]
%        [-n nav] [-bcmOuv] [-C coil_list] [-x wxoff] [-y wyoff] raw_file
%        isl1, isl2 are start and end slices to do
%        b0map = 0 for no shim correction
%              = 1 for normal correction [default]
%              = 2 for Bo map only
%        freq_offset (Hz) adds to b0map value
%        imgtyp is a bit map: 1=mag, 2=phase, 4=real, 8=imag
%        npix is power of two output matrix size
%        gamma is in Hz/G [default 4257] 
%        gain is factor for image scaling [default 1] 
%        zoomf is magnification factor
%        fermif is filter radius/kmax, 0 for none [0]
%        nav = 1 for navigator correction [(views==1)? 1:0]
%        -b to swap bytes in output images
%        -c to save coil images in multicoil
%        -m for no maxwell field correction
%        -O to generate images in single output file
%        -u to generate images for unfold
%        -v to print slice/echo notifications (verbose is slower)
%        coil_list is concatenated list of coils to use,
%            e.g., -C 0147 for coils 0, 1, 4 & 7
%        wxoff is offset for field map xslope
%        wyoff is offset for field map yslope
% This is source rev 85: 10/5/2008  signa rev 14.5

function [fmap,brain] = reconstructspiralfieldmap(filename_fMRI)

% function [fmap,brain] = reconstructspiralfieldmap(filename_fMRI)
%
% <filename_fMRI> is the location of the raw P-file with the spiral data
%
% return two volumes.  both volumes are double format and have dimensions 
% N x N x slices.
%  <fmap> is the fieldmap in Hz
%  <brain> is the second of the two volumes (magnitude brain)
%
% history:
% 2011/03/28 - comment out unnecessary section; don't write 'f' files; 
%              don't remove temp directory; don't write temporary dumb
%              files (at the expense of memory); code is much faster now!
% 2011/03/08 - Atsushi fixed the extra squaring bug.
%
% this routine requires the spirec executable.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% match the file
filename_fMRI = matchfiles(filename_fMRI);
assert(length(filename_fMRI)==1,'<filename_fMRI> does not match exactly one file');

% figure out a temporary directory
tempdir = tempname;
mkdirquiet(tempdir);

% get data to that directory
assert(copyfile(filename_fMRI{1},tempdir));

% temporarily change to the temp directory
olddir = pwd;
cd(tempdir);

% setup the variables
filename_fMRI = stripfile(filename_fMRI{1},1);  % just get the filename itself
filename_fieldmap = filename_fMRI;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% KK COMMENTED OUT:
% !ls -l P?????.7

spirec_options = ''; % '--multi 101';




% KK COMMENTED THIS OUT
% filename_fMRI     = input('Please enter the P-file name for the fMRI series : ','s');
% filename_fieldmap = input(sprintf('Please enter the P-file name for the field map [%s]  : ',filename_fMRI),'s');
% if (length(filename_fieldmap)==0);
%     filename_fieldmap = filename_fMRI;
% end
% 
% recon_res = input('What resolution do you want to reconstruct? [64,100,128,256, etc] ');
% if (size(recon_res) == size([]))
%    recon_res=128;
% end




%%%%%% Combine the field maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% KK COMMENTED OUT:
% map_res = recon_res; % 16; % recon_res;

% Clear the output files (had a problem when spirec overwrote links pointing to field map files)

disp(sprintf('!rm -f %s[._]*',filename_fMRI));
eval(sprintf('!rm -f %s[._]*',filename_fMRI));

% Get the in-plane off-isocenter locations from the header [THIS SECTION ADDED MAR 7, 2011]
[d,d,d,d,d,d,d,d,rhuser]=f_getspiralhdr(filename_fMRI);
dx = rhuser(23);
dy = rhuser(24);

  % KK REMOVED -m and -R from the following:
disp(sprintf('!spirec -u %g %g --fmaponly --savefmap -r %s -t %s',dx,dy,filename_fieldmap, filename_fMRI))
eval(sprintf('!spirec -u %g %g --fmaponly --savefmap -r %s -t %s',dx,dy,filename_fieldmap, filename_fMRI))

% M-file to read raw header and do sum of squares over coils for series of images
% filename_fieldmap = input('Please enter the P-file name for the series : ','s');
image_type  = 'freq'; % frequency maps
fid=fopen(filename_fMRI);
% Now check if the file exists
if (fid == -1) % did not find file so get required input
   nslc =    input('Please enter the number of slices :');
   necho =   input('                           echos  :');
   nfrms =   input('                           frames :');
   nrecv =   input('How many receivers did the scan have? ');
%%%if (nrecv == 1)
%%%    error('No need for sum of squares if you only have one receiver')
%%%end
else
   fclose(fid); % Found it so close it and read header
   [nslc,necho,nviews,nfrms,cols,rhptsize,rec_st,rec_end,rhuser]=f_getspiralhdr(filename_fieldmap);
%%%if (rec_st == rec_end),
%%%   error('We only had one receiver... exiting now');
%%%end
   if (necho > 1)
      error('I do not know how to handle more than one echo... exiting now');
   end
   nrecv = rec_end-rec_st+1;
   rows=nviews/rhuser(1)/2;
   k_limit=rhuser(17);
   g_max=rhuser(18)*10;              % mT/m
   s_max=rhuser(19);                 % SR mT/m/msec
   nturns=rhuser(20);
   gamma=rhuser(12);                 % Hz/G
   GAM=gamma*10.0;                   % to get mT/m
   tsp=rhuser(14)*1e-6;              % seconds
   interleaves=rhuser(1);
   sample_frequency=1/tsp;           % samples per second
   deltat=1/sample_frequency;
   t=(0:cols-1)*deltat;
end

% Let's figure out how large the images are...
% xyres = input('Please enter the matrix size (128, 256,etc) :');
% thisfilename = sprintf('%s.%s_%03dR%03dD%03d',filename_fMRI,image_type,slice,recnum,frame);
thisfilename = sprintf('%s.%s_%03d',filename_fMRI,image_type,0);
fid = fopen(thisfilename ,'r', 'ieee-le');
if (fid == -1)
   error(sprintf('Filename %s not found',thisfilename))
end
[thisfile, xyres2]=fread(fid,Inf,'float');
fclose(fid);

xyres = sqrt(xyres2); % assume a square matrix size

% loop over echos
% for echo = 0:necho-1
%   disp(sprintf('echo = %g',echo))
% loop over slices


% initialize (KK)
fmap = zeros(xyres,xyres,nslc);
brain = zeros(xyres,xyres,nslc);


   for slice =0:nslc-1
      disp(sprintf('slice = %g',slice))
% loop over frames to calculate mean
      frame = 0
      disp(sprintf('frame = %g',frame))
      disp('Calculating Mean....')     %%%% THIS LOOP IS SLOW, CAN WE SPEED UP?
      frames = 0;
      Sx   = zeros(xyres*xyres,1); % Sum of maps
      Sm   = zeros(xyres*xyres,1); % Sum of masks
% Loop over receivers
      for recnum=0:(nrecv-1)
% build up file name
         thisfilename = sprintf('%s.%s_%03d',filename_fMRI,image_type,recnum*nslc+slice);
         fid = fopen(thisfilename ,'r', 'ieee-le');
         if (fid == -1)
            error(sprintf('Filename %s not found',thisfilename))
         end
         thisfile=fread(fid,Inf,'float');
         fclose(fid);

         thismaskname = sprintf('%s.%s_%03d',filename_fMRI,'mask',recnum*nslc+slice);
         fid = fopen(thismaskname ,'r', 'ieee-le');
         if (fid == -1)
            error(sprintf('Filename %s not found',thismaskname))
         end
         thismask=fread(fid,Inf,'float');
         fclose(fid);
   
         frames = frames + 1;
%%%%%    Sx = Sx+thisfile.*thismask;
%%%%%    Sm = Sm + thismask;
         Sx = Sx+thisfile.*thismask;
         Sm = Sm + thismask;
      end
      mean = Sx./Sm;

% KK COMMENTED THIS OUT
%      imagesc(reshape(mean,xyres,xyres)');
%      colormap('gray')
%      title('Weighted Mean');, axis('square'), % pause(5)
%      drawnow, % pause(1)



% KK ADDED:
      fmap(:,:,slice+1) = reshape(mean,[xyres xyres]);
      brain(:,:,slice+1) = reshape(Sm.^(0.5),[xyres xyres]);


% KK COMMENTED THIS OUT, LET'S KEEP IT IN MEMORY ABOVE
% % Write out to file
%       outfilename = sprintf('%s.%s_%03d',filename_fMRI,'freq',slice);
% %     eval(sprintf('!mv %s %s.replaced',outfilename, outfilename))
%       fid=fopen(outfilename,'w','ieee-le');
%       if (fid == -1)
%          error(sprintf('Cannot open output file %s',outfilename))
%       end
%       fwrite(fid,mean,'float');
%       fclose(fid);
% %     disp(sprintf('Wrote to a file called %s',outfilename))
% 
%       outfilename = sprintf('%s.%s_%03d',filename_fMRI,'mask',slice);
% %     eval(sprintf('!mv %s %s.replaced',outfilename, outfilename))
%       fid=fopen(outfilename,'w','ieee-le');
%       if (fid == -1)
%          error(sprintf('Cannot open output file %s',outfilename))
%       end
%       fwrite(fid,Sm.^(0.5),'float');
%       fclose(fid);
% %     disp(sprintf('Wrote to a file called %s',outfilename))



% KK COMMENTED THIS OUT, WE DON'T ACTUALLY USE IT   
%       outfilename = sprintf('%s.%s_%03d',filename_fMRI,'f',slice);
% %     eval(sprintf('!mv %s %s.replaced',outfilename, outfilename))
%       fid=fopen(outfilename,'w','ieee-le');
%       if (fid == -1)
%          error(sprintf('Cannot open output file %s',outfilename))
%       end
%       fwrite(fid,floor(mean+0.5),'short');
%       fclose(fid);
% %     disp(sprintf('Wrote to a file called %s',outfilename))


%% KK COMMENTED THIS OUT, NOT NEEDED:
% % As a hack, create links to field map for each receiver
%       for recnum=1:(nrecv-1)
%          outfilename = sprintf('%s.%s_%03d',filename_fMRI,'freq',slice);
%          outlinkname = sprintf('%s.%s_%03d',filename_fMRI,'freq',recnum*nslc+slice);
% %         eval(sprintf('!ln -sf %s %s',outfilename, outlinkname))
%          cmd00 = sprintf('ln -sf %s %s',outfilename, outlinkname);
%          outfilename = sprintf('%s.%s_%03d',filename_fMRI,'mask',slice);
%          outlinkname = sprintf('%s.%s_%03d',filename_fMRI,'mask',recnum*nslc+slice);
% %         eval(sprintf('!ln -sf %s %s',outfilename, outlinkname))
%          cmd00 = [cmd00 '; ' sprintf('ln -sf %s %s',outfilename, outlinkname)];
%          outfilename = sprintf('%s.%s_%03d',filename_fMRI,'f',slice);
%          outlinkname = sprintf('%s.%s_%03d',filename_fMRI,'f',recnum*nslc+slice);
% %         eval(sprintf('!ln -sf %s %s',outfilename, outlinkname))
%          cmd00 = [cmd00 '; ' sprintf('ln -sf %s %s',outfilename, outlinkname)];
%          unix_wrapper(cmd00,0);  %% KK LUMPED IT ALL INTO ONE CALL
%       end % end of receiver loop


   end % end of slice
% end % end of echo





% KK COMMENTED THIS BIG SECTION OUT.
%
% %%%%% eval(sprintf('!spirec -l --multi 101 --loadfmap -r %s ',filename_fMRI))
% %%%%% eval(sprintf('!spirec -l --multi 101 --just 4 --loadfmap -r %s ',filename_fMRI))
% %%%%% eval(sprintf('!spirec -l --just 4 --loadfmap -r %s ',filename_fMRI))
% %%%%% eval(sprintf('!spirec -l --loadfmap -r %s -m 128 --just 2',filename_fMRI))
% %%%%% eval(sprintf('!spirec -l --loadfmap -r %s -m 128 --just 2 --multi 101',filename_fMRI))
% eval(sprintf('!spirec -l %s --loadfmap -r %s -m %d',spirec_options,filename_fMRI,recon_res))
% 
% 
% 
% % Do Sum of Squares on fMRI series
% 
% % M-file to read raw header and do sum of squares over coils for series of images
% %%%%% filename_fMRI = input('Please enter the P-file name for the series : ','s');
% %%%image_type  = input('Please enter the image type [I,Q,M or P]    : ','s');
% image_type  = 'M';
% fid=fopen(filename_fMRI);
% % Now check if the file exists
% if (fid == -1) % did not find file so get required input
%    nslc =    input('Please enter the number of slices :');
%    necho =   input('                           echos  :');
%    nfrms =   input('                           frames :');
%    nrecv =   input('How many receivers did the scan have? ');
% %  if (nrecv == 1)
% %      error('No need for sum of squares if you only have one receiver')
% %  end
% else
%    fclose(fid); % Found it so close it and read header
%    [nslc,necho,nviews,nfrms,cols,rhptsize,rec_st,rec_end,rhuser]=f_getspiralhdr(filename_fMRI);
% %  if (rec_st == rec_end),
% %     error('We only had one receiver... exiting now');
% %  end
%    nrecv = rec_end-rec_st+1;
%    rows=nviews/rhuser(1)/2;
%    k_limit=rhuser(17);
%    g_max=rhuser(18)*10;              % mT/m
%    s_max=rhuser(19);                 % SR mT/m/msec
%    nturns=rhuser(20);
%    gamma=rhuser(12);                 % Hz/G
%    GAM=gamma*10.0;                   % to get mT/m
%    tsp=rhuser(14)*1e-6;              % seconds
%    interleaves=rhuser(1);
%    sample_frequency=1/tsp;           % samples per second
%    deltat=1/sample_frequency;
%    t=(0:cols-1)*deltat;
% end
% 
% %  nfrms = 2;
% 
% % Let's figure out how large the images are...
% % xyres = input('Please enter the matrix size (128, 256,etc) :');
% thisfilename = sprintf('%s.%s_%03dR%03dD%03d',filename_fMRI,image_type,0,0,0);
% if (necho > 1)
%    thisfilename = sprintf('%sE%03d',thisfilename,echo);
% end
% 
% fid = fopen(thisfilename ,'r', 'ieee-le');
% if (fid == -1)
%    error(sprintf('Filename %s not found',thisfilename))
% end
% [thisfile, xyresfMRI2]=fread(fid,Inf,'short');
% fclose(fid);
% 
% xyresfMRI = sqrt(xyresfMRI2); % assume a square matrix size
% 
% if (xyres ~= xyresfMRI)
%     error('matrix size of fMRI does not match field map')
% end
% 
% % loop over echos
% for ech = 0:necho-1     % KK RENAMED echo to ech to avoid conflicts
%    disp(sprintf('echo = %g',ech))
% % loop over slices
%    for slice =0:nslc-1
%       disp(sprintf('slice = %g',slice))
% % loop over frames to calculate mean
%       for frame = 0:nfrms-1
%          disp(sprintf('frame = %g',frame))
%          disp('Calculation Sum of Squares....')
%          frames = 0;
%          Sx   = zeros(xyres*xyres,1);
% % Loop over receivers
%          for recnum=0:(nrecv-1)
% % build up file name
%             thisfilename = sprintf('%s.%s_%03dR%03dD%03d',filename_fMRI,image_type,slice,recnum,frame);
%             if (necho > 1)
%                thisfilename = sprintf('%sE%03d',thisfilename,ech);
%             end
%             fid = fopen(thisfilename ,'r', 'ieee-le');
%             if (fid == -1)
%                error(sprintf('Filename %s not found',thisfilename))
%             end
%             thisfile=fread(fid,Inf,'short');
%             fclose(fid);
% 
% % KK COMMENTED THIS OUT.  [WE GOING TO CLEAN UP OUR TEMP DIRECTORY ANYWAY]
% % SLOW
% % % Delete the individual receiver file image
% %             eval(sprintf('!rm %s',thisfilename))
%    
%             frames = frames + 1;
%             Sx = Sx+abs(thisfile).^2; % Take abs in case we load complex data
%          end
%          mean = sqrt(Sx)/frames;
% 
%          %% KK COMMENTED THIS OUT
%          %%imagesc(reshape(mean,xyres,xyres)');
%          %%colormap('gray')
%          %%title('Mean');, axis('square'), drawnow, % pause(1)
% 
% % Write out to file
%          outfilename = sprintf('%s.%s_%03dD%03d',filename_fMRI,image_type,slice,frame);
%          if (necho > 1)
%             outfilename = sprintf('%sE%03d',outfilename,ech);
%          end
%          fid=fopen(outfilename,'w','ieee-le');
%          if (fid == -1)
%             error(sprintf('Cannot open output file %s',outfilename))
%          end
%          fwrite(fid,mean,'short');
%          fclose(fid);
%          disp(sprintf('Writing to a file called %s',outfilename))
% 
%       end % end of frame
%    end % end of slice
% end % end of echo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KK ADDITIONS:

% KK CHANGED THINGS.  THIS IS NOW OBSOLETE:
% % load in results
% fmap = [];
% brain = [];
% for s=1:nslc
%   d = 0;
%   temp = double(loadbinary(sprintf('%s.freq_%03d',filename_fMRI,s-1),'float'));
%   temp2 = double(loadbinary(sprintf('%s.mask_%03d',filename_fMRI,s-1),'float'));
%   n = sqrt(length(temp)); assert(isint(n));
%   if s==1
%     fmap = placematrix2(zeros(n,n,nslc),reshape(temp,[n n]));
%     brain = placematrix2(zeros(n,n,nslc),reshape(temp2,[n n]));
%   else
%     fmap(:,:,s) = reshape(temp,[n n]);
%     brain(:,:,s) = reshape(temp2,[n n]);
%   end
% end

% change back to original directory and delete the temporary directory [ACTUALLY NO]
cd(olddir);
assert(rmdir(tempdir,'s'));





% JUNK:
%
%% <recon_res> is the resolution to reconstruct at
%%%%%%  brain(:,:,s) = double(loadbinary(sprintf('%s.mask_%03dD%03d',filename_fMRI,s-1,d),'int16',[##recon_res recon_res]));



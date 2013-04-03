function [timeframes,timekeys,digitrecord,trialoffsets] = ...
  ptviewmovie(images,frameorder,framecolor,frameduration,fixationorder,fixationcolor,fixationsize, ...
              grayval,detectinput,wantcheck,offset,moviemask,movieflip,scfactor,allowforceglitch, ...
              triggerfun,framefiles,frameskip,triggerkey,specialcon,trialtask)

% function [timeframes,timekeys,digitrecord,trialoffsets] = ...
%   ptviewmovie(images,frameorder,framecolor,frameduration,fixationorder,fixationcolor,fixationsize, ...
%               grayval,detectinput,wantcheck,offset,moviemask,movieflip,scfactor,allowforceglitch, ...
%               triggerfun,framefiles,frameskip,triggerkey,specialcon,trialtask)
% 
% <images> is a .mat file with 'images' as a uint8
%   A x B x 1/3 x N matrix with different images along the fourth dimension.  
%   the third dimension can be of size 1 or 3.  can be a cell vector, in 
%   which case we just concatenate all images together (the A and B 
%   dimensions should be consistent across cases).  images are referred 
%   to by their index (1-indexed).  if combining grayscale and color
%   images, we just repmat the grayscale images along the third dimension.
%   <images> can also be the uint8 matrix directly.  <images> can
%   also be a cell vector of uint8 elements, each of which is A x B x N 
%   or A x B x 3 x N (we detect this case if the size along the fourth
%   dimension is more than 1).  in these cases, the A and B must be 
%   consistent across the elements and <moviemask> must be [].
% <frameorder> (optional) is a vector of non-negative integers.  positive
%   integers refer to specific images in <images>.  zero means to show
%   no image (just the background and fixation).  this vector determines 
%   the order of presentation of images.  default is 1:N where N is the 
%   total number of images.
% <framecolor> (optional) is length(<frameorder>) x 3 with values in [0,255].
%   each row indicates a multiplication on color channels to apply for a 
%   given frame.  default is 255*ones(length(<frameorder>),3) which means 
%   to multiply each channel by 1 (i.e. do nothing special).  note that
%   when an entry in <frameorder> is 0, <framecolor> has no effect for that entry.
% <frameduration> (optional) is how many monitor refreshes you want a
%   single movie frame to last.  default: 15.
% <fixationorder> (optional) is:
%   (1) a vector of length 1+length(<frameorder>)+1 with elements that are non-
%       negative numbers in [0,1].  elements indicate alpha values.  the first 
%       element indicates what to use before the movie starts; the last element 
%       indicates what to use after the movie ends.  this is the regular
%       case, in which there is a single fixation dot color whose alpha value gets
%       modulated on different frames.
%   (2) elements can be negative integers plus an additional element that is an 
%       alpha value (thus, the length is 1+length(<frameorder>)+1+1).  the 
%       negative integers indicate row indices into <fixationcolor>, and the 
%       alpha value is used for all of the fixation colors.  in this alternative 
%       case, there are multiple possible fixation dot colors, all
%       of which get blended with a fixed alpha value.
%   (3) {A B C D E F} where A is the font size in (0,1), which is relative to
%       the size indicated by the first element of <fixationsize>; B is [ON OFF] 
%       where ON is a non-negative integer indicating the number of
%       frames for a digit to last and OFF is a non-negative integer
%       indicating the number of frames in between successive digits; C is whether
%       to omit the gray disc background; D is 0 (meaning just randomly pick), 1
%       (meaning randomly pick but ensure that successive digits are unique), or
%       a fraction between 0 and 1 indicating what proportion (on average) of cases
%       should have repeated digits (if you supply a negative fraction, we will
%       ensure that at most 2 successive digits will be the same); E 
%       is a positive integer indicating the number of CLUT 
%       entries to allocate at the end of the gamma table for pure
%       white and/or black (you need only one when F is 0, but you need two when
%       F is 1, and you can allocate more if you want to take up
%       entries); and F (optional) is 0 means white digits and 1 means alternate
%       white and black digits (default: 0).  we show a stream of random digits.
%       if C is 0, then the digits are shown on a gray disc that is
%       smoothly alpha-blended into the rest of the stimulus, with the
%       characteristics of this disc being determined by <fixationsize>. if C is
%       1, then the digits are directly superimposed on the rest of the stimulus.
%       before and after the movie, we show the digit '0'. note that
%       <fixationcolor> is ignored when <fixationorder> is of the {A B C D E F}
%       case.  also, note that if C is 1, then only the first element of
%       <fixationsize> is used (to determine the transparent box within
%       which the digits are presented).
%   default: ones(1,1+length(<frameorder>)+1).
% <fixationcolor> (optional) is a uint8 vector of dimensions 1 x 3, 
%   indicating the color to use for the fixation dot.  when <fixationorder>
%   is the special negative-integers case, <fixationcolor> should be a uint8
%   matrix of dimensions max(-fixationorder) x 3, indicating the possible colors
%   to use for the fixation dot.
%   default: uint8([255 255 255]).
% <fixationsize> (optional) is the size in pixels for the fixation dot.
%   default: 5.  special case is [A B] where A is the size for the fixation dot
%   and B is the number of pixels for the width of the border of the dot 
%   (the border is black).  for example, [6 1] means to use a dot that has
%   radius 3 and that has a border from 2 pixels to 3 pixels away from the center.
%   note that the default of 5 is equivalent to [5 0].
% <grayval> (optional) is the background color as uint8 1x1 or 1x3.
%   default: uint8(127).
% <detectinput> (optional) is whether to attempt to detect input during the 
%   showing of the movie.  if set to 1, you risk inaccuracies in 
%   the recorded times (for sure) and reduction (maybe) of your 
%   ability to achieve the desired framerate.
%   default: 1.
% <wantcheck> (optional) is whether to show some posthoc diagnostic figures
%   via ptviewmoviecheck.m.  default: 1.
% <offset> (optional) is [X Y] where X and Y are the
%   horizontal and vertical offsets to apply.  for example,
%   [5 -10] means shift 5 pixels to right, shift 10 pixels up.
%   default: [0 0].
% <moviemask> (optional) is an A x B matrix with values in [0,1].
%   0 means to pass through; 1 means to block.  we 
%   apply this mask to the images specified by <images>.  we blend 
%   with <grayval>.  default is [] which means do not apply a mask.
%   special case is when A or B (or both) are equal to 1; in this case,
%   we automatically expand (via bsxfun.m) to match the size of the images.
%   be careful: applying the mask is potentially a slow operation!
% <movieflip> (optional) is [J K] where J is whether to flip first
%   dimension and K is whether to flip second dimension.  this flipping
%   gets applied also to the fixation-related items.
%   default to [0 0].
% <scfactor> (optional) is a positive number with the scaling to apply
%   to the images in <images>.  if supplied, we multiply the number
%   of pixels in each dimension by <scfactor> and then round.  we use
%   bilinear filtering when rendering the images.  default is 1, and in
%   this case, we use nearest neighbor filtering (which is presumably
%   faster than bilinear filtering).
% <allowforceglitch> (optional) is
%   0 means do nothing special
%   [1 D] means allow keyboard input 'p' to force a glitch of duration D secs.
%     note that this has an effect only if <detectinput> is on.
%     forcing glitches is useful for testing purposes.
%   default: 0.
% <triggerfun> (optional) is the function to call right before we start the movie.
%   default is [] which means do not call any function.
%   if supplied, we create a 'trigger' event in <timekeys>, recording
%   the time of completion.
% <framefiles> (optional) is an sprintf string like '~/frame%05d.png'.  if supplied,
%   we write images containing the actual final frames shown on the display to 
%   the filenames specified by <framefiles>.  the files are 1-indexed, from 1 
%   through length(<frameorder>).  since writing to disk takes time, you may need 
%   to artificially increase <frameduration> to avoid glitches.  special case is
%   {A B} where A is like '~/frame%05d.png' and B is [R C] with the image dimensions
%   to crop to (relative to the center).  we make the parent directory for you
%   automatically (if necessary).  default: [].
% <frameskip> (optional) is a positive integer indicating how many frames to skip
%   when showing the movie.  for example, <frameskip>==2 means to show the 1st, 3rd,
%   5th, ... frames.  default: 1.  can also be 1/N for some positive integer N.
% <triggerkey> (optional) is
%   [] means any key can start the movie
%   X means if the first character of KbName(keyCode) where keyCode is obtained from
%     KbWait is X, then start the movie.  for example, you could pass in X as '5'.
%   default: [].
% <specialcon> (optional) is {A B C D}
%   where A is a Psychtoolbox calibration string, e.g. 'cni_lcd'
%         B is a vector of contrast values in (0,100], where the entries of this
%           vector matches the entries in <images>.  we determine the unique entries 
%           in this vector and do some precomputations based on that.
%         C is a N x 3 set of CLUT entries to use at the end of the gamma table,
%           e.g. for the fixation dot.  these are the linear values (before gamma correction).
%           we allocate 256-N entries to deal with the normal stimulus display.
%           thus, values in <images> should range from 0 through 255-N.  note that
%           when <fixationorder> is the {A B C D E F} case, we ignore the C input and always
%           use E CLUT entries at the end of the gamma table.
%         D is how many movie frames before a gamma change to attempt to do the 
%           gamma change.  the reason for this is that the gamma changes seem to take
%           a relatively long time and trying to do it at the last minute produces
%           weird glitching behavior.  note that because we have to do gamma changes
%           ahead of time, no gamma changes will be performed for the first D movie frames,
%           so don't expect any!
%   if supplied, we do special handling of the gamma table.  for example, 50 contrast
%   value means to restrict the range of the values in the gamma table to 0.25 through 
%   0.75.  we use a contrast value of 100 before and after the movie.  note that we load
%   the gamma table whenever a contrast change is needed.  we do not touch the gamma table 
%   when blank frames are shown.  it is unknown what value to use for D in order to avoid 
%   glitching behavior, so please test your movie.  (of course, the value you use for D 
%   should be compatible with your stimulus paradigm!)  if [] or not supplied, do nothing special.
% <trialtask> (optional) is {A B C D E F G H} where
%     A indicates the trial design, a matrix of dimensions T x F.  T corresponds to different trials.
%       F corresponds to the total number of frames in the movie.  each row should be 0 except 
%       for a consecutive string of 1s indicating the duration of the trial.  
%       trials should be mutually exclusive with respect to frames.
%     B is the fraction of trials (on average) that should present a dot.  we randomly flip a coin
%       for each trial to decide whether a dot is presented on that trial.
%     C is a cell vector of elements. each element should be a 2 x V matrix, each column 
%       indicating a valid location for the dot.  the units should be signed x- and y-coordinates 
%       in pixel units and is to be interpreted relative to the fixation location.  note that
%       these locations (and the <trialtask> stuff) are not affected by <movieflip>.
%     D is a vector 1 x T with the mapping from trials to the elements of C.
%     E is a uint8 vector of dimensions 1 x 3, indicating the color to use for the dot.
%     F is the size in pixels for the dot.
%     G is a non-negative number in [0,1] indicating the alpha value to use for the dot.
%     H is the positive number of frames to show the dot for.  should be less than or equal to
%       the shortest trial in A.
%   if supplied, we randomly choose trials to present a dot and during those trials we 
%   present the dot at a random location and at a random point in time.
%   if [] or not supplied, do nothing special.
%
% return <timeframes> as a 1 x length(<frameorder>) vector with the time of each frame showing.
%   (time is relative to the time of the first frame.)
% return <timekeys> as a record of the input detected.  the format is {(time button;)*}.
%   where button is a string (single button depressed) or a cell vector of strings (multiple
%   buttons depressed).  for regular button presses, the recorded time is what KbCheck returns.
%   the very first entry in <timekeys> is special ('absolutetimefor0') and indicates the absolute
%   time corresponding to the time of the first frame; all other time entries are relative to the
%   time of the first frame.
% return <digitrecord> as a 1 x length(<frameorder>) vector with the digit (0-9) shown on each
%   frame.  only the onsets of digits are recorded; the rest of the entries are NaN.  this input 
%   will be returned as [] if <fixationorder> is not the {A B C D E F} case.
% return <trialoffsets> as 2 x length(<frameorder>) with NaNs in all columns except those columns
%   corresponding to the presentation of the dot; for these columns, the first and second rows
%   have the x- and y-offsets of the dot in pixels, respectively.  note that in this case,
%   positive means to the right and to the top.
%
% Here are some additional notes:
% - What happens in the presentation of the movie:
%     First, we fill the background with gray and draw the fixation.
%     Then, we wait for a key from any keyboard (see <triggerkey>).
%       (You can toggle a safe mode by pressing '='.  In the safe mode,
%       nothing will happen until '=' is pressed again.)
%     Next, we wait until next vertical retrace and then issue 
%       <triggerfun> and proceed to show the movie.
%     In the movie, each frame either results in filling with gray (i.e.
%       when <frameorder> is 0) or results in showing an image,
%       and then the fixation is drawn.
%     Finally, we fill the background with gray and draw the fixation.
% - Before starting the movie presentation, we call Priority(MaxPriority) and hide the cursor.
%   At the end of the presentation, we restore the original priority and show the cursor.
% - We attempt to achieve frame-accurate showing of the movie.  If we glitch (i.e. we
%   show a frame too late), we attempt to catch up using a simple strategy --- try to show
%   the next frame at the ideal/perfect time.  However, there is no guarantee that
%   catching up will actually succeed (e.g. if the frame rate is really high).
%   Testing your movie is key.
% - In processing the images for the movie, we perform the following in order (if applicable):
%   Mask the image, flip the image, scale the dimensions of the image, and offset the image.
%   Masking involves converting the image to double, performing the mask, and converting
%   back to uint8, so beware of numerical precision issues.
% - We do not attempt to pre-call any functions.  You should make sure to perform a dry run of
%   your movie to make sure all mex files, functions, etc. are cached (for best performance).
% - During the movie presentation, we KbCheck all devices.  Note that KbCheck reports 
%   buttons for as long as they are held down.  The ESCAPE key forces an immediate exit 
%   from the movie.  
% - Even if <detectinput>, it is possible that there is no time to actually read input.
%   So it is important to test your particular setup!
%
% history:
% 2012/09/11 - add <trialtask>
% 2011/11/02 - add F to <fixationorder>; allow fourth argument of <fixationorder> to be negative
% 2011/11/02 - fix bug (would have crashed)
% 2011/10/26 - make <fixationorder> more flexible; make sure flipping happens for fixation stuff
% 2011/10/23 - tweak the <fixationorder> case to become {A B C D E}
% 2011/10/22 - add <fixationorder> {A B C D} case; add output <digitrecord>
% 2011/10/13 - add <specialcon>
% 2011/09/16 - add <triggerkey> and a = safe mode
% 2011/07/30 - add special entry for 'absolutetimefor0'
% 2011/04/22 - make parent directory of <framefiles> automatically.
% 2011/04/22 - fix <framefiles> bug by making sure the images are written out AFTER the flip.
%              (it previously was not writing out the last frame.)
% 
% example:
% pton;
% [timeframes,timekeys,digitrecord,trialoffsets] = ptviewmovie(uint8(255*rand(100,100,3,100)),[],[],2);
% ptoff;

% to do:
% - Rush? [probably not]
% - test: VBLSyncTest
% - useful for reference: [times,badout,misses]=CheckFrameTiming(100,1,10,5,0,0,0);
% - speedup: [resident [texidresident]] = Screen('PreloadTextures', windowPtr [, texids]);
% - run only in mirror mdoe? [probably not]
% - async flips was causing problems
% - using "don't clear" in flip was causing problems.

% input
if ischar(images)
  images = {images};
end
if ~exist('frameorder','var') || isempty(frameorder)
  frameorder = [];  % deal with later
end
if ~exist('framecolor','var') || isempty(framecolor)
  framecolor = [];  % deal with later
end
if ~exist('frameduration','var') || isempty(frameduration)
  frameduration = 15;
end
if ~exist('fixationorder','var') || isempty(fixationorder)
  fixationorder = [];  % deal with later
end
if ~exist('fixationcolor','var') || isempty(fixationcolor)
  fixationcolor = uint8([255 255 255]);
end
if ~exist('fixationsize','var') || isempty(fixationsize)
  fixationsize = 5;
end
if ~exist('grayval','var') || isempty(grayval)
  grayval = uint8(127);
end
if ~exist('detectinput','var') || isempty(detectinput)
  detectinput = 1;
end
if ~exist('wantcheck','var') || isempty(wantcheck)
  wantcheck = 1;
end
if ~exist('offset','var') || isempty(offset)
  offset = [0 0];
end
if ~exist('moviemask','var') || isempty(moviemask)
  moviemask = [];
end
if ~exist('movieflip','var') || isempty(movieflip)
  movieflip = [0 0];
end
if ~exist('scfactor','var') || isempty(scfactor)
  scfactor = 1;
end
if ~exist('allowforceglitch','var') || isempty(allowforceglitch)
  allowforceglitch = 0;
end
if ~exist('triggerfun','var') || isempty(triggerfun)
  triggerfun = [];
end
if ~exist('framefiles','var') || isempty(framefiles)
  framefiles = [];
end
if ~exist('frameskip','var') || isempty(frameskip)
  frameskip = 1;
end
if ~exist('triggerkey','var') || isempty(triggerkey)
  triggerkey = [];
end
if ~exist('specialcon','var') || isempty(specialcon)
  specialcon = [];
end
if ~exist('trialtask','var') || isempty(trialtask)
  trialtask = [];
end
if length(fixationsize)==1
  fixationsize = [fixationsize 0];
end
wantframefiles = ~isempty(framefiles);
if ~isempty(framefiles)
  if ischar(framefiles)
    framefiles = {framefiles []};
  end
end
if iscell(fixationorder) && (length(fixationorder) == 5)
  fixationorder{6} = [];
end
if iscell(fixationorder) && isempty(fixationorder{6})
  fixationorder{6} = 0;
end
if iscell(fixationorder) && ~isempty(specialcon)
  if fixationorder{6}==0
    specialcon{3} = repmat([255 255 255],[fixationorder{5} 1]);
  else
    specialcon{3} = repmat([255 255 255],[fixationorder{5} 1]);
    specialcon{3}(end-1,:) = [0 0 0];
  end
end

%%%%%%%%%%%%%%%%% DEAL WITH THE IMAGES

% load in the images
fprintf('loading images: starting...\n');
if iscell(images) && ischar(images{1})
  moviefile = images;
  images = [];
  for p=1:length(moviefile)
    temp = load(moviefile{p},'images');
    if ~isempty(images) & size(images,3)==1 & size(temp.images,3)==3  % do we need to make images into color?
      images = repmat(images,[1 1 3]);
    end
    if size(images,3)==3 & size(temp.images,3)==1  % do we need to make the temp.images into color?
      temp.images = repmat(temp.images,[1 1 3]);
    end
    images = cat(4,images,temp.images);
  end
  clear temp;
end
assert(isa(images,'uint8') || all(cellfun(@(x) isa(x,'uint8'),images)));  % check sanity
fprintf('loading images: done\n');

% deal with mask
if ~isempty(moviemask)
  fprintf('applying mask: starting...\n');
% OLD WAY:
%   moviemask = repmat(moviemask,[1 1 size(images,3) size(images,4)]);
%   images = uint8(moviemask * double(grayval) + (1 - moviemask) .* double(images));
%   clear moviemask;
  chunks = chunking(1:size(images,4),10);
  for p=1:length(chunks)  % to save on memory
    images(:,:,:,chunks{p}) = uint8(bsxfun(@plus,moviemask * double(grayval),bsxfun(@times,1 - moviemask,double(images(:,:,:,chunks{p})))));
  end
  fprintf('applying mask: done\n');
end

% deal with movieflip
if movieflip(1) && movieflip(2)
  wantflipfun = 1;
  flipfun = @(x) flipdim(flipdim(x,1),2);
elseif movieflip(1)
  wantflipfun = 1;
  flipfun = @(x) flipdim(x,1);
elseif movieflip(2)
  wantflipfun = 1;
  flipfun = @(x) flipdim(x,2);
else
  wantflipfun = 0;
end

%%%%%%%%%%%%%%%%% PREP

% get information about the PT setup
win = firstel(Screen('Windows'));
rect = Screen('Rect',win);

% calc
if iscell(images)
  d1images = size(images{1},1);
  d2images = size(images{1},2);
  if size(images{1},4) > 1
    dimwithim = 4;
  else
    dimwithim = 3;
  end
  csimages = cumsum(cellfun(@(x) size(x,dimwithim),images));
  nimages = csimages(end);
else
  d1images = size(images,1);
  d2images = size(images,2);
  nimages = size(images,4);
end

% deal with input (finally)
if isempty(frameorder)
  frameorder = 1:nimages;
end
if isempty(framecolor)
  framecolor = 255*ones(length(frameorder),3);
end
if isempty(fixationorder)
  fixationorder = ones(1,1+length(frameorder)+1);
end

% prepare movierect and fixationrect
movierect = CenterRect([0 0 round(scfactor*d2images) round(scfactor*d1images)],rect) + [offset(1) offset(2) offset(1) offset(2)];
fixationrect = CenterRect([0 0 2*fixationsize(1) 2*fixationsize(1)],rect) + [offset(1) offset(2) offset(1) offset(2)];  % allow doubling of fixationsize for room for anti-aliasing

% prepare fixation image
if ~iscell(fixationorder)

  fixationcase = any(fixationorder < 0);  % 0 means regular case, 1 means negative-integer case
    % 2*fixationsize x 2*fixationsize x 3 x N; several different uint8 solid colors
  fixationimage = zeros([2*fixationsize(1) 2*fixationsize(1) 3 size(fixationcolor,1)]);
  temp = find(makecircleimage(2*fixationsize(1),fixationsize(1)/2-fixationsize(2)));  % this tells us where to insert color
  for p=1:size(fixationcolor,1)
    temp0 = zeros([2*fixationsize(1)*2*fixationsize(1) 3]);  % everything is initially black
    temp0(temp,:) = repmat(fixationcolor(p,:),[length(temp) 1]);  % insert color in the innermost circle
    fixationimage(:,:,:,p) = reshape(temp0,[2*fixationsize(1) 2*fixationsize(1) 3]);
  %OLD:    fixationimage(:,:,:,p) = repmat(reshape(fixationcolor(p,:),[1 1 3]),[2*fixationsize 2*fixationsize]);
  end
  fixationalpha = 255*makecircleimage(2*fixationsize(1),fixationsize(1)/2);  % 2*fixationsize x 2*fixationsize; double [0,255] alpha values (255 in circle, 0 outside)

else

  % prepare digits as 2*fixationsize x 2*fixationsize x 3 x N; uint8 format
  digits = drawtexts(2*fixationsize(1),0,0,'Helvetica',fixationorder{1}, ...
                     [1 1 1],[0 0 0],mat2cell('0':'9',1,ones(1,10)));
  digits = round(normalizerange(digits,0,1));  % binarize so that values are either 0 or 1
  digsize = sizefull(digits,3);
  digits = repmat(vflatten(digits),[1 3]);  % T x 3
  if length(grayval)==1
    grayval0 = repmat(grayval,[1 3]);
  else
    grayval0 = grayval;
  end
  whzero = digits(:,1)==0;
  digits(whzero,:) = repmat(grayval0,[sum(whzero) 1]);  % 0 maps to the grayval
  digits(~whzero,:) = 255;  % 1 maps to the white value (remember to reserve this in the <specialcon> case)
  fixationimage = uint8(permute(reshape(digits,digsize(1),digsize(2),digsize(3),3),[1 2 4 3]));
    % make copy with black
  fixationimageB = fixationimage;
  fixationimageB(fixationimageB==255) = choose(isempty(specialcon),0,254);  % 1 maps to the black value (remember to reserve this in the <specialcon> case)
  fixationimage = cat(4,fixationimage,fixationimageB);
    % finally, add pure gray frame
  fixationimage(:,:,:,end+1) = repmat(reshape(grayval0,[1 1 3]),[size(fixationimage,1) size(fixationimage,2)]);
  
  % prepare alpha as 2*fixationsize x 2*fixationsize x 21; uint8 [0,255] alpha values
  if fixationorder{3}
    % the digits themselves are the 255 alpha values
    fixationalpha = uint8(reshape(255*double(~whzero),digsize(1),digsize(2),digsize(3)));
    fixationalpha = cat(3,fixationalpha,fixationalpha);
    fixationalpha(:,:,end+1) = 0;  % the last frame is completely transparent
  else
    % (255 in circle, 0 outside). gradual ramp.
    fixationalpha = repmat(uint8(255*makecircleimage(2*fixationsize(1), ...
                    fixationsize(1)/2-fixationsize(2),[],[],fixationsize(1)/2)),[1 1 21]);
  end
  
end

% prepare trial image
if ~isempty(trialtask)
  trialimage = repmat(reshape(trialtask{5},1,1,[]),[2*trialtask{6} 2*trialtask{6}]);
  trialalpha = uint8(trialtask{7} * (255*makecircleimage(2*trialtask{6},trialtask{6}/2)));
end

% now deal with flipping of fixation stuff
fixationimage = flipdims(fixationimage,movieflip);
fixationalpha = flipdims(fixationalpha,movieflip);

% prepare window for alpha blending
Screen('BlendFunction',win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

% init variables, routines, constants
timeframes = repmat(NaN,[1 floor((length(frameorder)-1)/frameskip)+1]);
timekeys = {};
digitrecord = [];
trialoffsets = [];
digitframe = [];
when = 0;
oldPriority = Priority(MaxPriority(win));
HideCursor;
mfi = Screen('GetFlipInterval',win);  % re-use what was found upon initialization!
filtermode = choose(scfactor==1,0,1);
getoutearly = 0;
glitchcnt = 0;
sound(0,1);

% precomputations for case when images is a cell of uint8
iscellimages = iscell(images);
if iscellimages
  whs = zeros(length(frameorder),2);
  for frame=1:length(frameorder)
    if frameorder(frame) ~= 0
      whs(frame,1) = firstel(find(frameorder(frame) <= csimages));
      whs(frame,2) = size(images{whs(frame,1)},dimwithim) - (csimages(whs(frame,1))-frameorder(frame));
    end
  end
end

% make directory
if wantframefiles
  mkdirquiet(stripfile(framefiles{1}));
end

% precompute cluts
if ~isempty(specialcon)

  % load and prep the cal
  [cal,cals] = LoadCalFile(specialcon{1});
  cal = SetGammaMethod(cal,0);

  % how many entries do we dedicate to the display of the images?
  nn = 256-size(specialcon{3},1);
  
  % what is a list of all the contrast levels that we will be using?
  allcons = union([100],specialcon{2});  % include 100 since we will use that before and after movie

  % do it
  specialcluts = [];
  for p=1:length(allcons)
    mn = 0.5 - 0.5*(allcons(p)/100);
    mx = 0.5 + 0.5*(allcons(p)/100);
    linearValues = [ones(3,1)*linspace(mn,mx,nn) specialcon{3}'];
    specialcluts(:,:,p) = PrimaryToSettings(cal,linearValues)';
  end

end

% deal with figuring out digit sequence for special fixation task
if iscell(fixationorder)

    % this will record onsets (for user consumption)
  digitrecord = NaN*zeros(1,ceil(length(frameorder)/sum(fixationorder{2})) * sum(fixationorder{2}));
    % this will tell us what to put on each frame [entries 1-10 map to '0':'9' white.
    %                                              entries 11-20 map to '0':'9' black.
    %                                              entry 21 maps to gray.]
  digitframe = digitrecord;
  digitpolarity = digitrecord;  % 0 means white.  1 means black.
  lastdigit = NaN;
  p = 1; cnt = 0; repeated = 0;
  while 1
    if p > length(frameorder)
      break;
    end
    if fixationorder{4}==1
      while 1
        digit = floor(rand*10);
        if ~isequal(digit,lastdigit)
          break;
        end
      end
    elseif fixationorder{4}==0
      digit = floor(rand*10);
    else
      case1 = fixationorder{4} < 0 && repeated;  % if we want to avoid triplets and we have just repeated...
      case2 = fixationorder{4} < 0 && ~repeated;  % if we want to avoid triplets and we haven't just repeated...
      case3 = fixationorder{4} > 0;  % if we just want to plow ahead
      if (case2 || case3) && (rand <= abs(fixationorder{4}))
        digit = lastdigit;
        repeated = 1;
        if isnan(digit)
          digit = 0;
          repeated = 0;
        end
      else
        repeated = 0;
        tttt = setdiff(0:9,lastdigit);
        ssss = randperm(length(tttt));
        digit = tttt(ssss(1));
      end
    end
    digitrecord(p) = digit;
    digitframe(p-1+(1:sum(fixationorder{2}))) = [repmat(digit+1,[1 fixationorder{2}(1)]) repmat(21,[1 fixationorder{2}(2)])];
    digitpolarity(p-1+(1:sum(fixationorder{2}))) = mod(cnt,2);
    lastdigit = digit;
    p = p + sum(fixationorder{2});
    cnt = cnt + 1;
  end
  
  % now truncate
  digitrecord = digitrecord(1:length(frameorder));
  digitframe = digitframe(1:length(frameorder));
  digitpolarity = digitpolarity(1:length(frameorder));
  
end

% figure out trialtask stuff
if ~isempty(trialtask)
  numtrials = size(trialtask{1},1);
  numframes = size(trialtask{1},2);
  dotrials = find(rand(1,numtrials) <= trialtask{2});  % indices of trials to show a dot on
  trialoffsets = NaN*zeros(2,numframes);  % compute x- and y-offsets for each frame
  for pp=1:length(dotrials)
    dotframes = find(trialtask{1}(dotrials(pp),:));  % indices of frames to show the dot on
      % choose a random duration within the trial by ignoring a random number of frames at the beginning
    dotframes = dotframes(floor(rand*(length(dotframes)-trialtask{8}+1)) + (1:trialtask{8}));
    locs = trialtask{3}{trialtask{4}(dotrials(pp))};  % 2 x L matrix of potential locations
    whichloc = ceil(rand*size(locs,2));  % pick one location. this is the index.
    trialoffsets(:,dotframes) = repmat(locs(:,whichloc),[1 length(dotframes)]);
  end
end

%%%%%%%%%%%%%%%%% START THE EXPERIMENT

% draw the background and fixation
Screen('FillRect',win,grayval,rect);
if iscell(fixationorder)
  texture = Screen('MakeTexture',win,cat(3,fixationimage(:,:,:,1),fixationalpha(:,:,1)));
else
  if fixationcase==0
    texture = Screen('MakeTexture',win,cat(3,fixationimage,uint8(fixationorder(1)*fixationalpha)));
  else
    texture = Screen('MakeTexture',win,cat(3,fixationimage(:,:,:,-fixationorder(1)),uint8(fixationorder(end)*fixationalpha)));
  end
end
Screen('DrawTexture',win,texture,[],fixationrect,[],0);
Screen('Close',texture);
if ~isempty(specialcon)
  Screen('LoadNormalizedGammaTable',win,specialcluts(:,:,allcons==100),1);  % use loadOnNextFlip!
  lastsc = 100;
end
Screen('Flip',win);

% wait for a key press to start
fprintf('press a key to begin the movie. (make sure to turn off network, energy saver, spotlight, software updates! mirror mode on!)\n');
safemode = 0;
while 1
  [secs,keyCode,deltaSecs] = KbWait(-3,2);
  temp = KbName(keyCode);
  if isequal(temp(1),'=')
    if safemode
      safemode = 0;
      fprintf('SAFE MODE OFF (the scan can start now).\n');
    else
      safemode = 1;
      fprintf('SAFE MODE ON (the scan will not start).\n');
    end
  else
    if safemode
    else
      if isempty(triggerkey) || isequal(temp(1),triggerkey)
        break;
      end
    end
  end
end

  % IS THE PREVIOUS LINE (RELATED TO KBWAIT) RELIABLE?  SHOULD WE USE SOMETHING DIFFERENT, LIKE:??
  % % just wait for any press
  %    % KbWait is unreliable probably would need a device input as well
  %    % but this is not an option (Psychtoolbox 1.0.5)
  %    % pause;
  %    iwait = 0;
  %    while ~iwait
  %        [~, ~, c] = KbCheck(-1);
  %        if find(c) == KbName('t')
  %            iwait = 1;
  %        end
  %    end

% wait until next vertical retrace (to reduce run-to-run variability)
Screen('Flip',win);

% issue the trigger and record it
if ~isempty(triggerfun)
  feval(triggerfun);
  timekeys = [timekeys; {GetSecs 'trigger'}];
end

% show the movie
framecnt = 0;
for frame=1:frameskip:length(frameorder)+1
  framecnt = framecnt + 1;
  frame0 = floor(frame);

  % we have to wait until the last frame is done.  so this is how we hack that in.
  if frame0==length(frameorder)+1
    while 1
      if GetSecs >= when
        getoutearly = 1;
        break;
      end
    end
  end

  % get out early?
  if getoutearly
    break;
  end

  % if special 0 case, just fill with gray
  if frameorder(frame0) == 0
    Screen('FillRect',win,grayval,movierect);

  % otherwise, make a texture, draw it at a particular position
  else
    if iscellimages
      if dimwithim==4   % THIS IS VERY VERY UGLY
        if wantflipfun
          texture = Screen('MakeTexture',win,feval(flipfun,images{whs(frame0,1)}(:,:,:,whs(frame0,2))));
        else
          texture = Screen('MakeTexture',win,images{whs(frame0,1)}(:,:,:,whs(frame0,2)));
        end
      else
        if wantflipfun
          texture = Screen('MakeTexture',win,feval(flipfun,images{whs(frame0,1)}(:,:,whs(frame0,2))));
        else
          texture = Screen('MakeTexture',win,images{whs(frame0,1)}(:,:,whs(frame0,2)));
        end
      end
    else
      if wantflipfun
        texture = Screen('MakeTexture',win,feval(flipfun,images(:,:,:,frameorder(frame0))));
      else
        texture = Screen('MakeTexture',win,images(:,:,:,frameorder(frame0)));
      end
    end
    Screen('DrawTexture',win,texture,[],movierect,0,filtermode,1,framecolor(frame0,:));
    Screen('Close',texture);
  end
  
  % draw the fixation
  if iscell(fixationorder)
    if fixationorder{6}==1
      if digitframe(frame) == 21
        whtodo = 21;
      else
        whtodo = digitframe(frame) + 10*digitpolarity(frame);
      end
    else
      whtodo = digitframe(frame);
    end
    texture = Screen('MakeTexture',win,cat(3,fixationimage(:,:,:,whtodo),fixationalpha(:,:,whtodo)));
  else
    if fixationcase==0
      texture = Screen('MakeTexture',win,cat(3,fixationimage,uint8(fixationorder(1+frame0)*fixationalpha)));
    else
      texture = Screen('MakeTexture',win,cat(3,fixationimage(:,:,:,-fixationorder(1+frame0)),uint8(fixationorder(end)*fixationalpha)));
    end
  end
  Screen('DrawTexture',win,texture,[],fixationrect,0,0);
  Screen('Close',texture);
  
  % draw the trial task dot
  if ~isempty(trialtask)
    if ~all(isnan(trialoffsets(:,frame0)))
      trialrect = CenterRect([0 0 2*trialtask{6} 2*trialtask{6}],rect) + ...
                  [offset(1) offset(2) offset(1) offset(2)] + ...
                  repmat(trialoffsets(:,frame0)' .* [1 -1],[1 2]);
                  % multiply y-coordinate by -1 because in PT, positive means down
      texture = Screen('MakeTexture',win,cat(3,trialimage,trialalpha));
      Screen('DrawTexture',win,texture,[],trialrect,0,0);
      Screen('Close',texture);
    end
  end

  % give hint to PT that we're done drawing
  Screen('DrawingFinished',win);
  
  % read input until we have to do the flip
  while 1
  
    % load the gamma table (for a future frame)
    if ~isempty(specialcon)
      frameL = frame0 + specialcon{4};
      if frameL <= length(frameorder)
        if frameorder(frameL)==0  % if blank frame, who cares, don't change
        else
          con = specialcon{2}(frameorder(frameL));
          if lastsc ~= con
%            sound(sin(1:100),1);
            Screen('LoadNormalizedGammaTable',win,specialcluts(:,:,allcons==con));  % don't use loadOnNextFlip!
            lastsc = con;
          end
        end
      end
    end

    % if we are in the initial case OR if we have hit the when time, then display the frame
    if when == 0 | GetSecs >= when
  
      % issue the flip command and record the empirical time
      [VBLTimestamp,StimulusOnsetTime,FlipTimestamp,Missed,Beampos] = Screen('Flip',win,when);
%      sound(sin(1:2000),100);
      timeframes(framecnt) = VBLTimestamp;

      % if we missed, report it
      if Missed > 0 & when ~= 0
        glitchcnt = glitchcnt + 1;
        didglitch = 1;
      else
        didglitch = 0;
      end
      
      % get out of this loop
      break;
    
    % otherwise, try to read input
    else
      if detectinput
        [keyIsDown,secs,keyCode,deltaSecs] = KbCheck(-3);  % all devices
        if keyIsDown

          % get the name of the key and record it
          kn = KbName(keyCode);
          timekeys = [timekeys; {secs kn}];

          % check if ESCAPE was pressed
          if isequal(kn,'ESCAPE')
            fprintf('Escape key detected.  Exiting prematurely.\n');
            getoutearly = 1;
            break;
          end

          % force a glitch?
          if allowforceglitch(1) && isequal(kn,'p')
            WaitSecs(allowforceglitch(2));
          end

        end
      end
    end

  end

  % write to file if desired
  if wantframefiles
    if isempty(framefiles{2})
      imwrite(Screen('GetImage',win),sprintf(framefiles{1},framecnt));
    else
      imwrite(uint8(placematrix(zeros([framefiles{2} 3]),Screen('GetImage',win))),sprintf(framefiles{1},framecnt));
    end
  end

  % update when
  if didglitch
    % if there were glitches, proceed from our earlier when time.
    % set the when time to half a frame before the desired frame.
    % notice that the accuracy of the mfi is strongly assumed here.
    when = (when + mfi / 2) + mfi * frameduration - mfi / 2;
  else
    % if there were no glitches, just proceed from the last recorded time
    % and set the when time to half a frame before the desired time.
    % notice that the accuracy of the mfi is only weakly assumed here,
    % since we keep resetting to the empirical VBLTimestamp.
    when = VBLTimestamp + mfi * frameduration - mfi / 2;  % should we be less aggressive??
  end
  
end

% draw the background and fixation
Screen('FillRect',win,grayval,rect);
if iscell(fixationorder)
  texture = Screen('MakeTexture',win,cat(3,fixationimage(:,:,:,1),fixationalpha(:,:,1)));
else
  if fixationcase==0
    texture = Screen('MakeTexture',win,cat(3,fixationimage,uint8(fixationorder(end)*fixationalpha)));
  else
    texture = Screen('MakeTexture',win,cat(3,fixationimage(:,:,:,-fixationorder(end-1)),uint8(fixationorder(end)*fixationalpha)));
  end
end
Screen('DrawTexture',win,texture,[],fixationrect,[],0);
Screen('Close',texture);
if ~isempty(specialcon)
  Screen('LoadNormalizedGammaTable',win,specialcluts(:,:,allcons==100),1);  % use loadOnNextFlip!
end
Screen('Flip',win);

%%%%%%%%%%%%%%%%% CLEAN UP

% restore priority and cursor
Priority(oldPriority);
ShowCursor;

% adjust the times in timeframes and timekeys to be relative to the first time recorded.
% thus, time==0 corresponds to the showing of the first frame.
starttime = timeframes(1);
timeframes = timeframes - starttime;
if size(timekeys,1) > 0
  timekeys(:,1) = cellfun(@(x) x - starttime,timekeys(:,1),'UniformOutput',0);
end
timekeys = [{starttime 'absolutetimefor0'}; timekeys];

% report basic timing information to stdout
fprintf('we had %d glitches!\n',glitchcnt);
dur = (timeframes(end)-timeframes(1)) * (length(timeframes)/(length(timeframes)-1));
fprintf('projected total movie duration: %.10f\n',dur);
fprintf('frames per second: %.10f\n',length(timeframes)/dur);

% do some checks
if wantcheck
  ptviewmoviecheck(timeframes,timekeys);
end










% JUNK
%   fixation to draw.  0 means do not draw anything.  1 means use the color
%   in the first row of <fixationcolor>, 2 means use the color in the
%   second row of <fixationcolor>, and so forth. 
%  Screen('FillOval',win,fixationcolor(fixationorder(1),:),fixationrect);
%    Screen('FillOval',win,fixationcolor(fixationorder(1+frame),:),fixationrect);
%  Screen('FillOval',win,fixationcolor(fixationorder(end),:),fixationrect);

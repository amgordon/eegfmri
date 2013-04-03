function constructstimuli

% function constructstimuli
%
% this is not really a function per se, but rather, a set of examples
% that demonstrate the construction of various stimuli.  if you run 
% this function, we create lots of figure windows with various stimuli.
%
% beware of the rendering of gray levels (as opposed to pure black or pure white)
% --- you should always inspect your stimuli for quality after making them!

% define
res = 64;  % image resolution desired
n = 16;    % number of images to make (if applicable)

% UNIFORM WHITE NOISE
f = rand(res,res,n);
viewimages(f);

% BINARY WHITE NOISE
f = double(rand(res,res,n) > 0.5);
viewimages(f);

% GAUSSIAN WHITE NOISE
f = normalizerange(randn(res,res,n),0,1,-3,3,1,1);
viewimages(f);

% PINK NOISE
f = normalizerange(generatepinknoise(res,1,n,1),0,1,-3,3,1,1);
viewimages(f);

% NATURAL IMAGE PATCHES
% see preparemcgillimages.m

% NATURAL IMAGES
f = [];
for p=1:n
  imfile = sprintf('/research/stimuli/natrev/photos.redogray/image%06d.png',p);  % NEED TO STICK THIS SOMEWHERE GENERAL
  f(:,:,p) = (double(imread(imfile))/255) .^ 2;
end
if res ~= 500
  f = processmulti(@imresize,f,[res res],'lanczos3');
end
viewimages(f);

% OBJECTS
f = reshape(loadmulti('/research/stimuli/kriegeskorte/images.mat','images')',[175 175 92]) .^ 2;  % NEED TO STICK THIS SOMEWHERE GENERAL
f = processmulti(@imresize,f(:,:,1:n),[res res],'lanczos3');
viewimages(f);

% GRATINGS
f = makegratings2d(res,2.^(0:floor(log2(res/2))),8,4)/2 + .5;
viewimages(f);

% CHECKERBOARDS
f = drawcheckerboards(res,1/3,2,3,1,-4,1,1,[0 0 0],[1 1 1]);
viewimages(f);

% POLYGON
f = drawclosedcontours(res,0,0,.5,0,[0 0 0],0,[],[1 1 1],encapsulate(@coordpolygon,[3:8 100]));
viewimages(f);

% POLYGON OUTLINES
f = drawclosedcontours(res,0,0,.5,0,[],500*1/20,[0 0 0],[1 1 1],encapsulate(@coordpolygon,[3:8 100]));
viewimages(f);

% ORIENTED LINES
f = drawbars(res,1/10,Inf,1,0,8,[0 0 0],[1 1 1]);
viewimages(f);

% CIRCULAR GRATINGS
f = permute(mean(reshape(makegratings2d(res,2.^(0:floor(log2(res/2))),180,1)/2 + .5,res,res,180,[]),3),[1 2 4 3]);
viewimages(f);

% SECOND-ORDER GRATINGS
gr = makegratings2d(res,2.^(0:floor(log2(res/2))),8,4)/2 + .5;
wn = rand(res,res,size(gr,3));
f = wn .* gr + .5 * (1-gr);
viewimages(f);

% SQUARE GRATINGS
f = drawcheckerboards(res,1/2,2,1,0,8j,4,1,[1 1 1],[0 0 0]);
viewimages(f);

% SURROUND-SUPPRESSION GRATINGS
circ = drawclosedcontours(res,0,0,1/3,0,[0 0 0],0,[],[1 1 1],encapsulate(@coordpolygon,[360]));
grat = makegratings2d(res,res/8,8,1)/2 + .5;
f = repmat(grat(:,:,1) .* (1-circ),[1 1 8]) + grat .* repmat(circ,[1 1 8]);
viewimages(f);

% ORIENTED BARS
f = drawbars(res,1/5,2,3,3,8,[0 0 0],[1 1 1]);
viewimages(f);

% ANGLES
f = drawclosedcontours(res,0,0,.8,0,[],-500*1/20,[0 0 0],[1 1 1],encapsulate(@(x) coordangle(0,x),linspacecircular(0,2*pi,16)));
viewimages(f);

% LETTERS
f = drawtexts(res,0,-.015,'Helvetica',.5,[0 0 0],[1 1 1],mat2cell('A':'Z',1,ones(1,26)));
viewimages(f);

% BAR TEXTURE
f = [];
for d=[0 1 2 3]
  f = cat(3,f,drawbartexture(res,1/64,1/16,pi/8,17,d*1/128,[0 0 0],[1 1 1]));
end
viewimages(f);








% JUNK:
%
% LUMINANCE
% % luminance modulation
% stim{10} = repmat(reshape(colonalt3(0,1,100),[1 1 100]),[64 64 1]);
% 
% % make lum
% lum = repmat(reshape(colonalt3(0,1,3),[1 1 3]),[64 64 1]);
% 
% % write lum
% for p=1:3
%   imwrite(uint8(matrixnormalize(lum(:,:,p),0,255,0,1)),sprintf('LUM%d.png',p));
% end
% 
% 
% CONNOR STIMULI
% connor stimclass book chapter

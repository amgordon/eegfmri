function m = loadbinary(file,precision,msize,lastdimrange,dim)

% function m = loadbinary(file,precision,msize,lastdimrange,dim)
%
% <file> is a pattern matching one or more files (see matchfiles.m)
% <precision> is something like 'int16'
% <msize> (optional) is the expected dimensions of the matrix for one file.
%   one of the dimensions can be 0, in which case we figure out what
%   that number should be.  default: [1 0].
% <lastdimrange> (optional) is [A B] where A<=B and A and B are
%   indices referring to the last dimension of <m> (i.e. the
%   last entry in <msize>).  we read and return only this portion
%   of the matrix.  default is [], which means return the whole matrix.
% <dim> (optional) is the dimension along which to concatenate 
%   matrices from different files.  default: 1.
%
% read <file> and return a matrix.  we assume that each file contains
% some data (i.e. it isn't empty).  for the machine format, we assume 
% IEEE floating point with little-endian byte ordering ('l'; see fopen).
%
% see also savebinary.m.
%
% example:
% savebinary('test','uint8',repmat(0:255,[2 1]));
% isequal(loadbinary('test','uint8',[0 256],[255 256]),repmat([254 255],[2 1]))

% constants
machineformat = 'l';

% input
if ~exist('msize','var') || isempty(msize)
  msize = [1 0];
end
if ~exist('lastdimrange','var') || isempty(lastdimrange)
  lastdimrange = [];
end
if ~exist('dim','var') || isempty(dim)
  dim = 1;
end

% get file name
file = matchfiles(file);
assert(length(file) >= 1,'<file> does not match at least one file');

% loop through files
for p=1:length(file)

  % open file
  fid = fopen(file{p},'r',machineformat);
  assert(fid ~= -1,'<file> could not be opened for reading');
  
  % handle case of no specific range
  if isempty(lastdimrange)
  
    % read it all in
    m0 = fread(fid,Inf,['*' precision],0,machineformat);
  
    % figure out msize
    if any(msize==0)
      assert(sum(msize==0)==1);  % make sure only one is 0
      msize(msize==0) = prod(size(m0))/prod(msize(msize~=0));  % calc appropriate value
      assert(all(isint(msize)),'<msize> is not correct');
    end
  
    % reshape
    m0 = reshape(m0,msize);
  
  % handle case of specific range
  else
  
    % peek to find out the byte size of a word
    temp = fread(fid,1,['*' precision],0,machineformat);
    wordsize = getfield(whos('temp'),'bytes');
    
    % check how big the data segment is
    assert(fseek(fid,0,'bof')==0);
    pos1 = ftell(fid); assert(pos1~=-1);
    assert(fseek(fid,0,'eof')==0);
    pos2 = ftell(fid); assert(pos2~=-1);
    
    % calculate number of words
    numwords = (pos2-pos1)/wordsize; assert(isint(numwords));
    
    % figure out msize
    if any(msize==0)
      assert(sum(msize==0)==1);  % make sure only one is 0
      msize(msize==0) = numwords/prod(msize(msize~=0));  % calc appropriate value
      assert(all(isint(msize)),'<msize> is not correct');
    end
    
    % calc slice size in terms of number of words
    slicesize = prod(msize(1:end-1));
    
    % calc number of slices wanted
    numslices = lastdimrange(2)-lastdimrange(1)+1;
    
    % read data and reshape
    assert(fseek(fid,slicesize*(lastdimrange(1)-1)*wordsize,'bof')==0);
    m0 = reshape(fread(fid,numslices*slicesize,['*' precision],0,machineformat),[msize(1:end-1) numslices]);
  
  end
  
  % close file
  assert(fclose(fid)==0);
  
  % save
  if p==1
    m = m0;  % get the datatype right instead of initializing to double via []
  else
    m = cat(dim,m,m0);
  end

end

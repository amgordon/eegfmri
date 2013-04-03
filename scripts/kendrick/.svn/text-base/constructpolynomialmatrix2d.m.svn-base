function [f,str] = constructpolynomialmatrix2d(matrixsize,locs,degree)

% function [f,str] = constructpolynomialmatrix2d(matrixsize,locs,degree)
%
% <matrixsize> is a 2D matrix size like [100 50]
% <locs> is a row or column vector of indices into that matrix size
% <degree> is the maximum polynomial degree desired
% 
% return <f>, a matrix of dimensions length(<locs>) x N
% with polynomial basis functions evaluated at <locs> in
% the columns.  the polynomial basis functions are evaluated
% over the range [-1,1] which is presumed to correspond to
% the beginning and ending element along each of the two dimensions.
% (if a dimension has only one element, the values are all set to 1.)
% also, return <str>, the algebraic expression that corresponds to
% the columns of <f>.  'x' refers to the first matrix dimension; 'y'
% refers to the second matrix dimension.
%
% note that there may be various gain factors on the basis functions
% (e.g. don't assume that they are all unit-length).
%
% also, beware of numerical precision issues for high degrees...
%
% see also constructpolynomialmatrix3d.m.
%
% example:
% [f,str] = constructpolynomialmatrix2d([30 30],find(ones(30,30)),2);
% str
% f = reshape(f,30,30,[]);
% figure; imagesc(makeimagestack(f));

% prep
x = sym('x');
y = sym('y');

% do the algebra
str = char(expand((x+y+1)^degree));

% sort the stuff in between + signs to try to ensure consistent ordering!!!
str = sort(strsplit(str,'+'));
str = cat(1,str,repmat({'+'},[1 length(str)]));
str = cat(2,str{:});
str = str(1:end-1);

% add a little padding so the 1 step below will work for degree 0
str = [' ' str ' '];

% change * to .*
  old = pwd;  % THIS IS A TOTAL HACK TO AVOID FUNCTION NAME CONFLICTS!
  cd(fullfile(matlabroot,'toolbox','matlab','funfun'));
str = vectorize(str);
  cd(old);

% remove +
str(str=='+') = ' ';

% change 1 to ones(size(x),1)
str0 = strrep(str,' 1 ',' ones(size(x)) '); assert(length(str0) ~= length(str));
str = str0;

% add brackets
str = [ '[' str ']' ];

% prep the linear coordinates
[x,y] = ind2sub(matrixsize,locs(:));
if matrixsize(1)~=1
  x = normalizerange(x,-1,1,1,matrixsize(1));
end
if matrixsize(2)~=1
  y = normalizerange(y,-1,1,1,matrixsize(2));
end

% do it
f = eval(str);

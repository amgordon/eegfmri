function f = nonlinearprfmodel(pp,dd,res,xx,yy,mode)

% function f = nonlinearprfmodel(pp,dd,res,xx,yy,mode)
%
% <pp> is the parameters as [R C SD GAIN N]
% <dd> is the stimuli as data x px*px
% <res> is pixels determining the resolution
% <xx>,<yy> (optional) are speed-ups.  pre-compute them like this:
%   [d,xx,yy] = makegaussian2d(res,2,2,2,2);
% <mode> (optional) is
%   0 means compute the response the fast and easy way.  however, this will choke
%     on the case of very small exponents (e.g. 0.003 or less) and small Gaussians,
%     in the sense that weird discontinuities may arise.
%   1 means to compute the responses with special handling to avoid the artifact
%     that may arise in <mode>==0.  the strategy is to look only at the non-zero
%     stimulus values (in <dd>) and temporarily boost the values in the Gaussian
%     up to exp(500).  by doing so, we help avoid the zero problems.  unfortunately,
%     this method is quite slow (seems to be about 20 times slower than the original
%     method).
%   default: 0.
%
% the output of this function is the response of the PRF model
% as a column vector of dimensions data x 1.
%
% note that small std dev for the Gaussian risks NaNs in the output.
% this is because the Gaussian may evaluate to all zeros.
%
% history:
% 2011/07/16 - first official version

% internal notes:
% - there is no nanreplace in this function.  this is because fitprf and fitprfstatic
%   automatically do the nanreplace check for us.
% - danger zones:
%   - if mode is 0, then really compressive exponents (e.g. 0.003 or less) combined
%     with small Gaussians may cause discontinuities in responses...
% - TODO: can we speed up <mode>==1?

% internal constants
goto = 500;  % what is the maximum exp exponent value to attain (in the mode==1 case)?

% input
if ~exist('xx','var') || isempty(xx)
  [d,xx,yy] = makegaussian2d(res,2,2,2,2);
end
if ~exist('mode','var') || isempty(mode)
  mode = 0;
end

% do it
switch mode
case 0

  % ah, this is an easy one-liner!
  f = pp(4) * ((dd*vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),pp(3),xx,yy,0,0))).^pp(5));

case 99  % THIS MODE IS FOR INTERNAL USE ONLY!

  % ah, this is an easy one-liner!
  f = pp(4) * divnorm((dd*vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),pp(3),xx,yy,0,0))) / (2*pi*pp(3)*pp(3)),pp(5),pp(6));

case 999  % THIS MODE IS FOR INTERNAL USE ONLY!

  % ah, this is an easy one-liner!
  f = pp(4) * divnormpp((dd*vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),pp(3),xx,yy,0,0))) / (2*pi*pp(3)*pp(3)),pp(5),pp(6),pp(7));

case 1

  % do the gaussian but omit the exponentiation.  this is row vector of exponent values (e.g. -100).
  gau = flatten(makegaussian2d(res,pp(1),pp(2),pp(3),pp(3),xx,yy,0,1));
  
  % compute the dot-product with a boost to avoid zero problems
  %   notes: precaching the ~=0 and dd(:,good) does not speed up that much (about 3s to 2s)
  f =     zeros(size(dd,1),1);
  boost = zeros(size(dd,1),1);
  for p=1:size(dd,1)
    good = dd(p,:)~=0;      % logical vector indicating non-zero stimulus values [MOST EXECUTION TIME]
    temp = gau(good);       % what are the gaussian exponents we have to worry about?
    if isempty(temp)
      continue;  % if the stimulus is all zeros, leave f and boost set to 0
    end
    mx = max(temp);         % what is the maximum exponent that we have to worry about?
    boost(p) = goto-mx;     % how much exponent is it okay to boost by?
    f(p) = dd(p,good) * exp(temp + boost(p)).';   % do the boost, then do the dot-product [MOST EXECUTION TIME]
  end
      % EXPERIMENTAL
      % %   gaut = repmat(gau,[size(dd,1) 1]);
      % %   gaut(dd==0) = 0;
      % %   boost = max(gaut,[],2);
      % %   newgau = exp(repmat(boost,[1 size(gaut,2)])+gaut);
      % % %  newgau(isinf(newgau)) = 0;
      % %   f = dot(dd,newgau,2);
  
  % ok, do the PRF exponent, and then we have to remove the effect of the boost
  f = f.^pp(5) ./ exp(boost*pp(5));
  
  % finally, we have to apply the PRF gain parameter
  f = pp(4) * f;
 
end



%%%%%%%%%%%%%%%%%%%%%%%% JUNK:

%OLD  f = (pp(4) * fct.^pp(5)) * f;
%   if mode==1
%     % we want the Gaussian to sum to 1.  so we must apply this scaling to the Gaussian.
%     % note that if the Gaussian is all zeros (e.g. because the std dev is too small), then the final
%     % output will be non-finite.  this is okay.
% %OLD    fct = 1/l1vectorlength(exp(gau));
% %    gau = gau - log(l1vectorlength(exp(gau)));
%   else
% %OLD    fct = 1;
%   end
%  , as well as the special sum-to-one factor we computed earlier.
  
% - we now normalize the Gaussian to sum to one.  this makes the interpretation of
%   the gain parameter clearer (i.e. it's the response to a full-field stimulus).
%   - if sd of the Gaussian gets too small, the normalization of the Gaussian
%     will fail (divide by 0, leading to Inf).
%   - also, if the Gaussian roams far away, that also risks failure of normalization.
%   - these failures are okay, because the nanreplace check will ensure that if the
%     optimizer encounters these cases that they will show up as really bad fits and
%     so the optimizer will avoid them.


% JUNK:
%   % precompute
%   if ~iscell(dd)
%   
%     % first element (dd{1}) is the stimulus itself
%     dd = {dd};
%     
%     % second element (dd{2}) is logical indicating where non-zero values are
%     dd{2} = dd{1}~=0;
%     
%     % third element (dd{3}) is a cell vector of the non-zero values (as a row vector)
%     dd{3} = {};
%     for p=1:size(dd{1},1)
%       dd{3}{p} = dd{1}(p,dd{2}(p,:));
%     end
%   
%   end

% A function that can compute the various gradients used by classifierLogisticRegressionSimple,
% used inside the optimization loop in that function
%
% History:
% - 2009 June - created from previous code - fpereira@princeton.edu 
%

function [f,df,fpart,lambdapart] = classifierLogisticRegressionSimple_gradients(Wvectorized,X,Y,lambda,yindices,ymask,regularization,optimization);

[n,m] = size(X); k = size(ymask,2); tm = m + 1;
W = reshape(Wvectorized,[tm k]);

df = zeros(tm,k);

% 1) the E part, each example contributes to df estimates for its own class only

df(1,:) = sum(ymask,1);
for j = 1:k
  df(2:end,j) = sum(X(yindices{j},:),1)';
end

% 2) the log part (hence the - signs)

tmp1 = X * W(2:end,:) + repmat(W(1,:),n,1);
tmp2 = exp(tmp1); % n x k
tmp3 = tmp2 ./ repmat(sum(tmp2,2),1,k); % n x k

df(1,:) = df(1,:) - sum(tmp3,1);
for j = 1:k
  df(2:end,j) = df(2:end,j) - sum(X .* repmat(tmp3(:,j),[1,m]),1)';
end

%f_eterm   = sum(tmp1.*ymask,2);
%f_logterm = log(sum(tmp2,2));
%f  = sum( fterm - f_logterm);
f = sum( sum(tmp1.*ymask,2) - log(sum(tmp2,2)) );

fpart = f;

tmp = W(2:end,:);

% norms are not over bias term
switch regularization
 case {'L2'}
  lambdapart = 0.5*lambda*sum(tmp(:).^2); 
  f  = f - lambdapart;
  df(2:end,:) = df(2:end,:)-lambda*tmp; % vectorize it
end
df = df(:); % vectorize it


switch optimization
 case {'minimize'}
  % the function (-1 is because the optimization package does minimization)
  % (this is the only part of the bound where E appears, so only compute this)
  f  = -1*f;
  df = -1*df;
 otherwise
  % do nothing to the function and gradient output
end

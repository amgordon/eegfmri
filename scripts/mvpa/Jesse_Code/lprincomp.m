


function [coeff, score, latent] = lprincomp(x, max_dim)
% 
% function [coeff, score, latent] = lprincomp(x)
% 
% identical to matlab princomp but can handle very large matrices in an efficient way.
% 
% Inputs
%   x           a 2D matrix with rows corresponding to observations
%   max_dim     maximum number of principal components that should be calculated
%
% Outputs
%   read pincomp help for an the explanation of outputs
% 
% Developed by Roozbeh Kiani, May 27 2011
% 

if nargin<2
    dim = min(size(x));
else
    dim = min(max_dim, min(size(x)));
end

    %using pca function for large matrices
n = size(x,1);
        %center A
D = x-repmat(mean(x),[n 1]);
        %run svd
[U,S,V] = pca(D, dim);
coeff = V;
score = U.*repmat(diag(S)',[n 1]);
latent = diag(S).^2/(n-1);






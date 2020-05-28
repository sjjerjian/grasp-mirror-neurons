function [score,expl,alind,alinds,C] = subproj(X,V,k)
% [score,expl,alind,alinds,C] = subproj(X,V,k)
% Projects data X onto subspace defined by first k axes of V (principal axes)
% Assumes data is already centered !
%
% Returns projections of data onto first k axes of V
% Calculates explained variance of first k PCs
% Calculates alignment index (see below for explanation)
%
% If V is not provided or provided as empty, it will be default to the eigenvectors of cov(X).
% if k exceeds the number of axes (columns) provided in V , the code will just use all columns of V.
%
% INPUTS: 
% X     : data matrix (N*P matrix, observations x features)
% V     : projection axes (P*m matrix where m>=k)
% k     : number of axes to use (default is 5)
%
% OUTPUTS:
% score : projections of X onto k axes of V (N*length(k) matrix)
% expl  : explained variances of k PCs in subspace 1 and subspace2 (k x 2 matrix)
% alind : alignment index (scalar between 0 and 1)
% alinds: alignment along each axis (lower limit of 0, but can be >1 for individual axes)
% C     : covariance matrix of X
%
% Steven Jerjian 06/04/2018
%
%% Some checks
if nargin < 3, k = 5; end

if exist('V','var') && ~isempty(V) && k > size(V,2)
    disp('You asked for more PCs than there are dimensions of data...using all dimensions');
end

%% Project data into space, and calculate latent variables

C = cov(X); 
[v2,latent2,e] = pcacov(C);

% V defaults to the axes best aligned with input data
if nargin < 2 || isempty(V), V = v2; end 
k = min(k,size(V,2));

V       = V(:,1:k);    % cut down number of axes
score   = X * V;       % project the data onto these axes
latent1 = V' * C * V;  % latent variables in subspace defined by V

%% Calculate variance captured by each component 

% if data is centered this should be equivalent
% origVar   = sum(sum(bsxfun(@minus,X,mean(X)).^2));
% expl(:,1) = 100 * (sum(score.^2) ./ origVar);

expl(:,1) = 100 * (diag(latent1) ./ sum(latent2));  % variance explained by each PC relative to overall variance in data
expl(:,2) = e(1:k);                                 % variance explained in data

%% Alignment Index
%{
The alignment index is the ratio between the data variance explained in the chosen subspace to the
overall variance in the data (up to the first d dimensions)
if the subspace is defined by the eigenvectors of the data itself, then the alignment will be 1, and if
the data is fully orthogonal to the chosen subspace, then alignment will be 0.
%}

% trace(P1'*cov(D2)*P1) / trace(P2'*cov(D2)*P2) 
% equivalent to...
% trace(P1'*cov(D2)*P1) / eigs(cov(D2))

% where D2 = data, P1 are first k principal axes of one subspace, P2 are first k  principal axes of 2nd
% subspace, eigs are first k largest eigenvalues

alind = trace(latent1) / sum(latent2(1:k));

% also equivalent to...
% latent2 = v2' * C * v2;       % if V == v2, then latent1 = latent2, alind = 1
% alind = trace( latent1 ) / trace( latent2(1:k,1:k) );

% alignment for each axis independently (lower bound is zero, but no theoretical upper bound)
alinds = diag(latent1) ./ latent2(1:k); 

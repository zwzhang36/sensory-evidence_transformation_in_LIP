function [E,D] = pca_selfdefine(X)
%
% Inputs:
%   X  - [nb x nd] raw data
%           nb: number of observation
%           nd: number of dimension, 
% Outputs:1
%   E - [nd x nd] : eigenvector of the covariance matrix
%       E(1,:): the first component
%   D - [nd x nd] : eigenvaule


%% pca
%Normalization
X = X-ones(size(X,1),1)*mean(X);

% Calculate the eigenvalues and eigenvectors of the new covariance matrix.
covarianceMatrix = X'*X/size(X,1);
[E, D] = eig(covarianceMatrix); % eigenvectors as column

% Sort the eigenvalues  and recompute matrices
[~,order] = sort(diag(-D));
E = E(:,order);
E = E';
d = diag(D); 
D = diag(d(order));

% E = E';
% dsqrtinv = real(d.^(-0.5));
% Dsqrtinv = diag(dsqrtinv(order));
% V = Dsqrtinv*E;% V whitening matrix
% Y = E*X;%Y the project matrix of the input data X without whiting

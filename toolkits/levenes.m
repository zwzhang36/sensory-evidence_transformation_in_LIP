function [p,atab] = levenes(data,grp)
%
% function [p,atab] = levenes(data,grp)
%
% Implements Levene's test for the mogeneity of variance.
% Can be used on any data suitable for an n-way ANOVA.
%
% Parameters:
%   data and grp are identical to the first two parameters passed to
%   MatLab's anovan function. Type 'help anovan' for details.
%   
%   data should be a vector with measurement data (the dependent variable).
%   grp should be a 1 x nVar cell array, where nVar is the number of
%   variables (or factors). Each cell of grp should be a vector of same
%   size as data used as a grouping variable (as in anovan.m).
%
% Returns:
%   p and atab are the same as the first two values returned from anovan.
%   p is the probability of a difference in variance among the different
%   cells in the data matrix greater than or equal to the one in data.
%   atab is the usual information from an ANOVA, eg.:
%   
%    'Source'    'Sum Sq.'    'd.f.'    'Singular?'    'Mean Sq.'    'F'         'Prob>F'
%    'X1'        [ 0.0416]    [   3]    [        0]    [  0.0139]    [0.2797]    [0.8400]
%    'Error'     [30.3020]    [ 612]    [        0]    [  0.0495]          []          []
%    'Total'     [30.3436]    [ 615]    [        0]            []          []          []
%
% Notes:
%   The solution here is based on the fact that Levene's test is equivalent
%   to the following:
%      Let D = abs(y(i,j) - y_hat(.,j))
%      Where i is the data index and j indexes cell membership, and y_hat
%      is the within-cell mean ("cell" here refering to the cells of a data
%      matrix of).
%      
%      Then perform a one-way ANOVA with D as the dependent variable.
%
% Aaron Schurger
% June 2007
%

nVar = length(grp);
grpn = cell(size(grp));


% parse grouping variables into number sequence 
for i=1:nVar
    grpn{i} = zeros(size(grp{i}));
    u = unique(grp{i});
    for j = 1:length(u)
        ind = grp{i}==u(j);
        grpn{i}(ind) = j-1; % zero-based
    end
end

grpn_mat = cell2mat(grpn);
base = max(max(grpn_mat)+1); % bases, as in base-2, base-10, etc...
bases = base .^ [nVar-1:-1:0];

lgrp = (grpn_mat * bases') + 1; % back to 1-based
ind = unique(lgrp);
nGrp = length(ind); % number of cells (i.e. groups)

% now we're ready to calculate means for each group

M = zeros(size(data));
N = length(data);
for i=1:nGrp
    M(lgrp == ind(i)) = mean(data(lgrp == ind(i)));
end

D = abs(data - M); % deviation scores

[p,atab] = anovan(D,{lgrp},'model','full','display','off');

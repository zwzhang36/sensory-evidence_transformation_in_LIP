function [negLP,grad,H] = l1smoothness(prs, negloglifun, beta, num_reg)
% [negLP,grad,H] = neglogposterior(prs,negloglifun,Cinv)
%
% Compute negative log-posterior given a negative log-likelihood function
% and zero-mean Gaussian prior with inverse covariance 'Cinv'.
%
% Inputs:
 %   prs [d x 1] - parameter vector
%    negloglifun - handle for negative log-likelihood function
%   Cinv [d x d] - response (spike count per time bin)
%
% Outputs:
%          negLP - negative log posterior
%   grad [d x 1] - gradient 
%      H [d x d] - Hessian (second deriv matrix)

% Compute negative log-posterior by adding quadratic penalty to log-likelihood
prs_diff = diff(prs); % constant term 
cum_num  = cumsum([1 num_reg(1:end-1)]);
prs_diff(cum_num) = 0;

add_grad  = zeros(size(prs)); % gradient of each parameters
add_hess  = zeros(numel(prs), numel(prs)); % hessian matrix of each parameters
first_prs = cumsum([2 num_reg(1:end-1)]);
last_prs  = cumsum(num_reg)+1;
for i = 2:numel(prs)
    % the first parameter of a regressor group
    if any(i == first_prs)
        add_grad(i) = prs(i)*sign(prs(i)-prs(i+1));
        add_hess(i,i)   = sign(prs(i)-prs(i+1));
        add_hess(i,i+1) = 0;
        continue
    end
    % the last parameter of a regressor group
    if any(i == last_prs)
        add_grad(i) = prs(i)*sign(prs(i)-prs(i-1));
        add_hess(i,i)   = sign(prs(i)-prs(i-1));
        add_hess(i,i-1) = 0;
        continue
    end
    % the middle parameters of a regressor group
    add_grad(i) = prs(i)*sign(prs(i)-prs(i-1)) + prs(i)*sign(prs(i)-prs(i+1));
    add_hess(i,i)   = sign(prs(i)-prs(i-1)) + sign(prs(i)-prs(i+1));
    add_hess(i,i+1) = 0;
    add_hess(i,i-1) = 0; 
end

%%
switch nargout

    case 1  % evaluate function
        negLP = negloglifun(prs) + sum(abs(prs_diff))*beta;
    
    case 2  % evaluate function and gradient
        [negLP,grad] = negloglifun(prs);
        negLP = negLP + sum(abs(prs))*beta;        
        grad = grad + beta*add_grad;

    case 3  % evaluate function and gradient
        [negLP,grad,H] = negloglifun(prs);
        negLP = negLP + sum(abs(prs))*beta;        
        grad = grad + beta*add_grad;
        H = H + beta*add_hess;
end


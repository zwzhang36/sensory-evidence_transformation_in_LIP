function [negLP,grad,H] = l2Regl2Sms(prs, negloglifun, Cinv, beta, num_reg)
% Inputs:
%   prs [d x 1] - parameter vector
%   negloglifun - handle for negative log-likelihood function
%   Cinv [d x d] - response (spike count per time bin)
%
% Outputs:
%   negLP - negative log posterior
%   grad [d x 1] - gradient 
%   H [d x d] - Hessian (second deriv matrix)

% Compute negative log-posterior by adding quadratic penalty to log-likelihood

%%
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
        add_grad(i) = prs(i)-prs(i+1);
        add_hess(i,i)   = 1 - prs(i+1);
        add_hess(i,i+1) = prs(i) - 1;
        continue
    end
    % the last parameter of a regressor group
    if any(i == last_prs)
        add_grad(i) = prs(i)-prs(i-1);
        add_hess(i,i)   = 1 - prs(i-1);
        add_hess(i,i-1) = prs(i) - 1;
        continue
    end
    % the middle parameters of a regressor group
    add_grad(i) = 2*prs(i)-prs(i+1)-prs(i-1);
    add_hess(i,i)   = 2 - prs(i+1) - prs(i-1);
    add_hess(i,i+1) = 2*prs(i) - prs(i-1) - 1;
    add_hess(i,i-1) = 2*prs(i) - prs(i+1) - 1; 
end

%%
switch nargout

    case 1  % evaluate function
        negLP = negloglifun(prs) + .5*prs'*Cinv*prs + prs_diff'*prs_diff*beta;
    
    case 2  % evaluate function and gradient
        [negLP,grad] = negloglifun(prs);
        negLP = negLP + .5*prs'*Cinv*prs + prs_diff'*prs_diff*beta;        
        grad = grad + Cinv*prs + 2*beta*add_grad;

    case 3  % evaluate function and gradient
        [negLP,grad,H] = negloglifun(prs);
        negLP = negLP + .5*prs'*Cinv*prs + prs_diff'*prs_diff*beta;        
        grad = grad + 2*Cinv*prs + 2*beta*add_grad;
        H = H + Cinv + 2*beta*add_hess;
end


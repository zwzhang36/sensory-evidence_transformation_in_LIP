function mtest = modeltest(wml, varStd, X, Y, rglzt, smusns)
% % % get a summary for the regression
% rearrange the mean and std of the weight
% Inputs:
%   wml [d x 1] - parameter vector
%   varStd[d x 1] - handle for negative log-likelihood function
%   Y
%   X
%   rglzt:  the setting of the regularization
%   smusns: the setting of the smoothness

Y_hat = exp(X*wml);
SSresidual = mean((Y-Y_hat).^2); % mean squared error, GLM no offset
SStotal = mean((Y-mean(Y)).^2); % squared error of spike train

mtest.rsquare = 1-SSresidual/SStotal;
mtest.SSE = SSresidual;
mtest.SST = SStotal;

% log likelihood 
mtest.loglr = (Y'*log(Y_hat) -  sum(Y_hat));
mtest.loglr_avge  = mtest.loglr/numel(Y);

mtest.loglr_rglzt = mtest.loglr;
% penalty from regularization
if strcmp(rglzt.method, 'lasso')
    lambda = rglzt.lambda*numel(Y);
    mtest.loglr_rglzt = mtest.loglr_rglzt - sum(abs(wml))*lambda;
elseif strcmp(rglzt.method, 'ridge')
    lambda = rglzt.lambda*numel(Y);
    Cinv   = lambda*eye(numel(wml)); Cinv(1,1)=0;
    mtest.loglr_rglzt = mtest.loglr_rglzt - 0.5*wml'*Cinv*wml;
end

% penalty from smoothness
wml_diff = diff(wml); % constant term
cum_num  = cumsum([1 smusns.num_reg(1:end-1)]);
wml_diff(cum_num) = 0;
if strcmp(smusns.method, 'lasso')
    beta = smusns.beta*numel(Y);
    mtest.loglr_rglzt = mtest.loglr_rglzt - sum(abs(wml))*beta;
elseif strcmp(smusns.method, 'ridge')
    beta = smusns.beta*numel(Y);
    mtest.loglr_rglzt = mtest.loglr_rglzt - wml_diff'*wml_diff*beta;
end

% average prediction error
% mtest.deviation = norm(Y-pred_spike); for guessian distribution,not poisson

mtest.deviation = 2*(sum(Y.*log(Y./Y_hat + 1e-100)) - sum(Y) + sum(Y_hat));
mtest.deviation = mtest.deviation/numel(Y);

% baised estimation
mtest.pvalue = normcdf(0,abs(wml),varStd); % used for now, maybe not strict enough

end
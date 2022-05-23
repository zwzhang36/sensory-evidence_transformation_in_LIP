function [wml, varStd, nlogli] = ridgePoissonGLM(X,Y,lambda)

Imat = eye(size(X,2)); % identity matrix of size of filter + const
Imat(1,1) = 0; % remove penalty on constant dc offset
Cinv = lambda*Imat; % set inverse prior covariance

wInit = X \ Y;
fnlin = @linear_regression.expfun; % inverse link function (a.k.a. nonlinearity)
lossfunc = @(w)(linear_regression.poisson(w, X, Y, fnlin)); % cost/loss function
ridgelossfunc = @(w)linear_regression.negLogPosterior.l2regularization(w,lossfunc,Cinv);
% lfunc = @(w)(glms.neglog.bernoulli(w, XDegnMat, Y_spike)); % cost/loss function

opts = optimoptions(@fminunc, 'Display', 'off', 'Algorithm', 'trust-region', ...
    'GradObj', 'on', 'Hessian','on');
[wml, nlogli, exitflag, ostruct, grad, hessian] = fminunc(ridgelossfunc, wInit, opts);
% pGLMwts = glmfit(XDegnMat,Y_spike,'poisson', 'constant', 'off');
coVarMat = inv(hessian);% the invert hessian matrix is the covariance matrix
varStd = sqrt(diag(coVarMat));

end
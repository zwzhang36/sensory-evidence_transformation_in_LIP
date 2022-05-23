function [wml, varStd, nlogli] = penaltyPoissonGLM(X, Y, rgulzt, smusns, varargin)


num_reg = smusns.num_reg;

% estimate the initial weight
wInit = X \ Y;
fnlin = @linear_regression.expfun; % inverse link function (a.k.a. nonlinearity)
lossfunc = @(w)(linear_regression.poisson(w, X, Y, fnlin)); % cost/loss function

%% 
% for regularization
if isfield(rgulzt, 'lambda'); lambda = rgulzt.lambda*numel(Y); end
if strcmp(rgulzt.method ,'ridge')
    Imat = eye(size(X,2));    % identity matrix of size of filter + const
    Imat(1,1) = 0;            % remove penalty on constant dc offset
    Cinv_r = lambda*Imat;     % set inverse prior covariance
end
% for smoothness
if isfield(smusns, 'beta');   beta = smusns.beta*numel(Y);     end
if strcmp(smusns.method ,'ridge')
    Imat = eye(size(X,2));  % identity matrix of size of filter + const
    Imat(1,1) = 0;          % remove penalty on constant dc offset
    Cinv_s = beta*Imat;     % set inverse prior covariance
end


%% define the loss and gradient 
%  this part of code is not elegant, but faster 

if ~isempty(rgulzt.method) && isempty(smusns.method) 
    % only the penalize the magnitude of parameters
    if strcmp(rgulzt.method ,'lasso')
        penaltyloss = @(w)linear_regression.negLogPosterior.l1regularization(w,...
            lossfunc, lambda);
    elseif strcmp(rgulzt.method ,'ridge')
        penaltyloss = @(w)linear_regression.negLogPosterior.l2regularization(w,...
            lossfunc, Cinv_r);
    else
        error('unknown regularization method');
    end
elseif ~isempty(smusns.method) && isempty(rgulzt.method)
    % only the penalize the smoothness of parameters
    if strcmp(smusns.method ,'lasso')
        penaltyloss = @(w)linear_regression.negLogPosterior.l1smoothness(w,...
            lossfunc, beta, num_reg);
    elseif strcmp(smusns.method ,'ridge')
        penaltyloss = @(w)linear_regression.negLogPosterior.l2smoothness(w,...
            lossfunc, beta, num_reg);
    else
        error('unknown smoothness method');
    end
elseif ~isempty(rgulzt.method) && ~isempty(smusns.method)
    if strcmp(rgulzt.method, 'lasso') && strcmp(smusns.method, 'lasso')
        penaltyloss = @(w)linear_regression.negLogPosterior.l1Regl1Sms(w,...
            lossfunc, lambda, beta, num_reg);
    elseif strcmp(rgulzt.method, 'lasso') && strcmp(smusns.method, 'ridge')
        penaltyloss = @(w)linear_regression.negLogPosterior.l1Regl2Sms(w,...
            lossfunc, lambda, beta, num_reg);
    elseif strcmp(rgulzt.method, 'ridge') && strcmp(smusns.method, 'lasso')
        penaltyloss = @(w)linear_regression.negLogPosterior.l2Regl1Sms(w,...
            lossfunc, Cinv_r, beta, num_reg);
    elseif strcmp(rgulzt.method, 'ridge') && strcmp(smusns.method, 'ridge')
        penaltyloss = @(w)linear_regression.negLogPosterior.l2Regl2Sms(w,...
            lossfunc, Cinv_r, beta, num_reg);
    end
else 
    error('call the function without penalty');
end
penaltylossfunc = penaltyloss;
%%%%%%%%%%%%%%%%%
% add by zzw, 20210124
if ~isempty(varargin) && numel(varargin{:})~=0
    index = varargin{1, 1}.index+1;
    wmlGround = reshape([varargin{1, 1}.params.weight_mean],1,[]);
    varGround = reshape([varargin{1, 1}.params.weight_std],1,[]);
    wGround = reshape([varargin{1, 1}.params.weight_mean],1,[]);

    wInit(index) = wGround; % the first one is constant
    penaltylossfunc = @(w) penaltylossfuncGround(w, penaltyloss, index);
end
%%%%%%%%%%%%%%%%%

% lfunc = @(w)(glms.neglog.bernoulli(w, XDegnMat, Y_spike)); % cost/loss function
%% minimize the loss function
opts = optimoptions(@fminunc, 'Display', 'off', 'Algorithm', 'trust-region', ...
    'GradObj', 'on', 'Hessian','on');
[wml, nlogli, exitflag, ostruct, grad, hessian] = fminunc(penaltylossfunc,...
    wInit, opts);
% pGLMwts = lassoglm(X, Y, 'poisson','CV',10);

%%%%%%%%%%%%%%%%%
% add by zzw, 20210124
if ~isempty(varargin) && numel(varargin{:})~=0
    hessianNew = hessian;
    hessianNew(index,:) = [];hessianNew(:,index) = [];
    coVarMatNew = inv(hessianNew);% the invert hessian matrix is the covariance matrix
    varStdNew = sqrt(diag(coVarMatNew));
    
    numParams = size(hessian,1);
    varStd = zeros(numParams,1);
    indexNew = 1:numParams;
    for i = index
        indexNew(indexNew==i) = [];
    end
    varStd(index) = varGround;
    varStd(indexNew) = varStdNew;
%%%%%%%%%%%%%%%%%
else
    coVarMat = inv(hessian);% the invert hessian matrix is the covariance matrix
    varStd = sqrt(diag(coVarMat));
end
%{
% it is convenient for regularization (the parameter will not be necessary 
%  change with sampe size ) 

% Zhewei Zhang, 20210304, normalized by the sample size, 
% should be consistent with the function '+linear_regression.poisson.m'
% nomalize in both two function or neither.
varStd = varStd/numel(Y);
%}
end

function [L,dL,H] = penaltylossfuncGround(w, lossfunc, index)
    [L,dL,H] = lossfunc(w);
    dL(index) = 0; H(:,index) = 0;H(index,:) = 0;
end

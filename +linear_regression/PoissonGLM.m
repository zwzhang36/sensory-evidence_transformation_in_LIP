function [wml, varStd, nlogli] = PoissonGLM(X, Y, varargin)


wInit = X \ Y;
fnlin = @linear_regression.expfun; % inverse link function (a.k.a. nonlinearity)
loss = @(w)(linear_regression.poisson(w, X, Y, fnlin)); % cost/loss function
% lfunc = @(w)(glms.neglog.bernoulli(w, XDegnMat, Y_spike)); % cost/loss function
lossfunc = loss;
%%%%%%%%%%%%%%%%%
% add by zzw, 20210124
if ~isempty(varargin) && numel(varargin{:})~=0
    index = varargin{1, 1}.index+1;
    wmlGround = reshape([varargin{1, 1}.params.weight_mean],1,[]);
    varGround = reshape([varargin{1, 1}.params.weight_std],1,[]);
    wGround = reshape([varargin{1, 1}.params.weight_mean],1,[]);

    wInit(index) = wGround; % the first one is constant
    lossfunc = @(w) lossfuncGround(w, loss, index);
end
%%%%%%%%%%%%%%%%%
opts = optimoptions(@fminunc, 'Display', 'notify-detailed', 'Algorithm', 'trust-region', ...
    'GradObj', 'on', 'Hessian','on','OptimalityTolerance',1e-6);
    %,...'FunctionTolerance', 1e-6/numel(Y)
[wml, nlogli, exitflag, ostruct, grad, hessian] = fminunc(lossfunc, wInit, opts);
% pGLMwts = glmfit(XDegnMat,Y_spike,'poisson', 'constant', 'off');

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

end

function [L,dL,H] = lossfuncGround(w, lossfunc, index)
    [L,dL,H] = lossfunc(w);
    dL(index) = 0; H(:,index) = 0;H(index,:) = 0;
end

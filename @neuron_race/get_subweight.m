function theta = get_subweight(obj,type)
if nargin==2 && strcmp(type,'pre-assigned')
    theta = [0.9 0.54 0.3 -0.3 -0.54 -0.9 -0.9 -0.54 -0.3 0.3 0.54 0.9];
    obj.subweight = theta; return;
end

condition = cell2mat(obj.race.condition');
shapes = zeros(size(condition,1), 12);
for i = 1:12;  shapes(:,i) = sum(condition==i,2); end
X = shapes;
Y = obj.race.choice_color';
%%
if nargin==2 && strcmpi(type,'subjective_lasso')
    [B, fitinfo] = lassoglm(X,Y,'binomial','link','logit', ...
        'lambda', 2.7569e-04);
    theta = B;
    
    func = @(Y_hat) sum(Y.*log(Y_hat) + (1-Y).*(log(1-Y_hat)));
    Y_hat = exp(X*B + fitinfo.Intercept);
    Y_hat = Y_hat./(Y_hat+1);
    logL = func(Y_hat); numP = 13; numObs = size(X,1);
    [aicLasso,bicLasso] = aicbic(logL,numP,numObs,'Normalize',true);
    %{
                [B, fitinfo] = lassoglm(X,Y,'binomial','link','logit', ...
                    'CV',3);
                % lambda=2.1421e-05 for monkey M, population, MinDeviance
                % lambda=2.7569e-04 for monkey H, population, MinDeviance
                % lassoPlot(B, fitinfo,'plottype','CV'); legend('show')
                idxLambdaMinDeviance = fitinfo.IndexMinDeviance;
                theta = B(:,idxLambdaMinDeviance);
    %}
end

if nargin==2 && strcmp(type,'subjective')
    % logistic regression
    % mdl = fitglm(X,Y,'distribution','binomial');
    theta = glmfit(X,Y,'binomial','link','logit','constant','on');
    theta = theta(2:end);
end

if nargin==1 || (nargin==2 && strcmp(type,'symmetric'))
    % default setting
    % use the subjective weight from each day
    X = shapes(:,1:6) - shapes(:,7:12);
    % logistic regression
    mdl = fitglm(X,Y,'Distribution','binomial');
    % theta = glmfit(X,Y,'binomial','link','logit','constant','on');
    theta = mdl.Coefficients.Estimate;
    theta = [theta(2:end);-theta(2:end)];
    SE = mdl.Coefficients.SE;
    SE = [SE(2:end); SE(2:end)];
    pVal = mdl.Coefficients.pValue;
    pVal = [pVal(2:end); pVal(2:end)];
    
    % numP = 7; numObs = size(X,1);
    % [aicSymtc,bicSymtc] = aicbic(mdl.LogLikelihood, numP, numObs);
    %{
    X = shapes;
    % function used for calculating the logL in symmetric and slope model
    help = @(a,b,c) exp(X(:,1:6)*a + X(:,7:12)*a*b + c);
    pred = @(a,b,c) help(a,b,c)./(1+help(a,b,c));
    loss = @(p) -sum(  Y.*log(   pred(p(1:6),p(7),p(8))) + ...
                    (1-Y).*log(1-pred(p(1:6),p(7),p(8))));
    params = randn(8,1);
    theta = fmincon(loss, params);
    theta = [theta(1:6); -theta(1:6); -theta(7)];
    %}
end

%%
if diff(size(theta))<0; theta = theta'; end
obj.subweight = theta;
obj.subweight_SE = SE;
obj.subweight_pVal = pVal;
end
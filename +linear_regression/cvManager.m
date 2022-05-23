function [cvError, lambda_best, beta_best, results] = cvManager(nCV,...
    regularization, smoothness, designMat,sel_var, spCount, trials, num_reg)
import linear_regression.*
% assume the effect of regularization and smoothness are independent

%% input inspection

if nCV == 1
    error('what are you doing? No.CV==1?');
end

if isempty(regularization.lambda)
    regularization.lambda = 0;
%     warning('plz set the available parameters in the regularization');
end

if isempty(smoothness.beta)
    smoothness.beta = 0;
%     warning('plz set the available parameters for smoothness');
end

%% data set
nTrials = numel(trials);
suflInd = randperm(nTrials); % shuffle the trial order
degMat = designMat(sel_var);
degMat = [degMat{:}];

%% result saving

rglzt  = regularization;
smusns = smoothness;
if strcmpi(regularization.lambda,'search')
    lamdList = [0 10.^(-12:2)];
else
    lamdList = regularization.lambda;
end
if strcmpi(smoothness.beta,'search')
    betaList = [0 10.^(-12:2)];
else
    betaList = smoothness.beta;
end
nlamd = numel(lamdList);
nbeta = numel(betaList);

%%
results = cell(nlamd, nbeta); 
cvError = zeros(nlamd, nbeta)+nan;

if strcmpi(regularization.lambda,'search')
    nthLmda = 9; % initial value, 1e-6
else
    nthLmda = 1;
end

if strcmpi(regularization.lambda,'search')
    nthBeta = 11; % initial value, 1e-6
else
    nthBeta = 1;
end
% test
rglzt.lambda = lamdList(nthLmda); % initial value
smusns.beta  = betaList(nthBeta); % initial value

[result, divs] = glmCV(degMat, spCount, trials, suflInd, num_reg,...
    nCV, smusns, rglzt);

cvError(nthLmda, nthBeta) = mean(divs); % deviation
results{nthLmda, nthBeta} = result;
if numel(regularization.lambda) == 1 && numel(smoothness.beta)==1
    lambda_best = regularization.lambda;
    beta_best = smoothness.beta; 
    return 
end

while true
    params = [nthLmda,   nthBeta+1; nthLmda,   nthBeta-1;...
              nthLmda+1, nthBeta;   nthLmda-1, nthBeta];
    params(min(params')<1,:) = []; % not exceed the parameter list
    params(params(1,:)>nlamd, :) = [];
    params(params(2,:)>nbeta, :) = [];
    
    % test the soround parameters
    for param = params'
        if ~isnan(cvError(param(1), param(2))); continue; end
        
        rglzt.lambda = lamdList(param(1)); % initial value
        smusns.beta  = betaList(param(2)); % initial value        

        [result, divs] = glmCV(degMat, spCount, trials, suflInd, num_reg,...
            nCV, smusns, rglzt);
        
        cvError(param(1), param(2)) = mean(divs);
        results{param(1), param(2)} = result;
    end
    
    % stop if the cvError is smaller in both direction
    if cvError(nthLmda, nthBeta) == nanmin(cvError(:))
        break
    end
    
    % otherwise. use the minimum as now source
    [~, ndx] = nanmin(cvError(:));
    nthBeta = ceil(ndx/nlamd);         % nth smoothness.beta is the best
    nthLmda = ndx-(nthBeta-1)*nlamd; % nth regularization.lambda is the best
    
end

lambda_best = lamdList(nthLmda);
beta_best   = betaList(nthBeta);
fprintf('chosen lamdba:%0.2d || bata:%0.2d \n', lambda_best, beta_best);
end

function [results, cvError] = glmCV(degMat, spCount, trials, suffle, num_reg,...
    nCv, smusns, rglzt)

import linear_regression.*

cvError = zeros(nCv,1);
results = cell(nCv,1);
nTrials = numel(trials);

parfor nth = 1:nCv
    %% prepare data for training and test
    indStart = 1+round((nth-1)*nTrials/nCv);
    indStop  = round(nth*nTrials/nCv);
    
    test = zeros(nTrials,1)==1;
    test(indStart:indStop) = true;
    test = test(suffle);
    
    trianingset_degMat = cell2mat(degMat( trials(~test),:));
    trianingset_ySpike = cell2mat(spCount(trials(~test)))';
    
    testingset_degMat = cell2mat(degMat( trials(test),:));
    testingset_ySpike = cell2mat(spCount(trials(test)))';
    
    [trianingset_degMat, constindex, nocon_No] = ...
        preprocessing(trianingset_degMat,num_reg);
    [testingset_degMat,~,~] = preprocessing(testingset_degMat, num_reg);
    
    %% fitting and statistic test
    % training
    if rglzt.lambda==0 && smusns.beta==0
        [wml,varStd, nlogli] = PoissonGLM(trianingset_degMat,...
            trianingset_ySpike);
    else
        [wml,varStd, nlogli] = penaltyPoissonGLM(trianingset_degMat,...
            trianingset_ySpike, rglzt, smusns);
    end
    % testing
    result = modeltest(wml, varStd, testingset_degMat,...
        testingset_ySpike,  rglzt,  smusns);
    
    % save results
    weight_mean = repack(wml,    num_reg, constindex, nocon_No);
    weight_std  = repack(varStd, num_reg, constindex, nocon_No);
    result.weight_mean = weight_mean;
    result.weight_std  = weight_std;
    result.lambda = rglzt.lambda;
    result.beta = smusns.beta;
    
    cvError(nth) = result.deviation;
    results{nth} = result;
    
end

end
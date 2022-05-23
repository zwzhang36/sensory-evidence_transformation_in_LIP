function varargout = regression(obj,varargin)
% % % preform the regression
% % cross validtaion
% % regulation, ridge, l2 norm
import linear_regression.*
disp('GLM regression')
%% check the input arguement
if isempty(obj.designMat);error('design matrix is not caculated');end
[cv, regularization, smoothness, sel_var, conditionOn] = input_judge(varargin, obj.var, obj.order);
% trials for analysis
eval(['trials = ',obj.trials_equ,';']);
if int2unit(obj.neuron.spCount.race_fix_unit) ~= int2unit(obj.unit)
    obj.neuron.timepoint2count_fixation('interval',obj.unit);
    warning('the timebin for the timeline and spCount should be the same')
end
spikeCount = obj.neuron.spCount.race_fix.all;

num_reg = cell2mat(obj.var.num_regress(sel_var));

%% if validation ie needed, do all the regressions, return the result,
% compare the logLR to find the best parameters
if cv~=0
    [cvError, lambda_best, beta_best, result_cv] = cvManager(cv, regularization,...
        smoothness, obj.designMat, sel_var, spikeCount, trials, num_reg);
end

%% refit the model using all the data and best lambda
% seperate the data into training set and testing set
dataset = dataprepare(obj.designMat,sel_var,spikeCount,trials);
% preprocessing
[dataset.degMat, constindex, nocon_No] = preprocessing(dataset.degMat, num_reg);

if ~isempty(regularization.method)
    if cv==0; lambda_best = regularization.lambda; end
    lambdaList = regularization.lambda;
    regularization.lambda = lambda_best;
    regularization.lambda_best = lambda_best;
end

if ~isempty(smoothness.method)
    if cv==0; beta_best  = smoothness.beta; end
    betaList = smoothness.beta;
    smoothness.beta = beta_best;
    smoothness.beta_best = beta_best;
end

if isempty(regularization.method) && isempty(smoothness.method)
    % no penalty for coefficients
    [wml,varStd,~] = PoissonGLM(dataset.degMat, dataset.Y_spike, conditionOn);
else
    [wml,varStd,~] = penaltyPoissonGLM(dataset.degMat, dataset.Y_spike,...
        regularization, smoothness);    
end
result = modeltest(wml, varStd, dataset.degMat, dataset.Y_spike,...
    regularization, smoothness);
if ~isempty(smoothness.method);     smoothness.beta = betaList;         end
if ~isempty(regularization.method); regularization.lambda = lambdaList; end

% statitics test
weight_mean = repack(wml,    num_reg, constindex, nocon_No);
weight_std  = repack(varStd, num_reg, constindex, nocon_No);
result.pvalue = repack(result.pvalue, num_reg, constindex, nocon_No);
result.weight_mean = weight_mean;
result.weight_std  = weight_std;

if ~isempty(regularization.method)
    result.lambda = lambda_best;
end
if ~isempty(smoothness.method)
    result.beta   = beta_best;
end
%% summary the info about the neuron and variable

% info for neuron, trials and regressors
neuron = [obj.neuron.basic_info.monkey,'_',num2str(obj.neuron.basic_info.date),...
    '_',num2str(obj.neuron.basic_info.order)];
result.neuron = neuron;
result.trials = obj.trials_equ;
labels = fieldnames(obj.order);
result.label = {'constant',labels{sel_var}};
result.variable = obj.var;
result.smoothness = smoothness.method;
result.regularization = regularization.method;

result.cv = cv;
if cv~=0 
    result.result_cv = result_cv;
    result.cv_error = cvError;
end
prev_reg = length(obj.result);
obj.result{prev_reg+1} = result;
if nargout>=1
    varargout{1} = result;
end
if nargout>=2
    varargout{2} = result_cv;
end
end


function [cv,regularization,smoothness,sel_var, conditionOn] = input_judge(arginput,var,order)
% % % check the input
% % % cv: whether cross validation is needed
% % % regularization: ridge or lasso, for now, ridge only
% % % sel_var: the variable used for further analysis
cv = 0;
sel_var = [];
conditionOn = [];
regularization = struct('method','','lambda','');
smoothness = struct('method','','beta','');

%%
cv_label = find(strcmpi(arginput,'cv')==1)+1;
if ~isempty(cv_label)
    cv = arginput{cv_label};
    arginput(cv_label-1:cv_label) = [];
end

regu_label = find(strcmpi(arginput,'regularization')==1)+1;
if ~isempty(regu_label)
    regularization.method = arginput{regu_label};
    arginput(regu_label-1:regu_label) = [];
end

lambda_label = find(strcmpi(arginput,'lambda')==1)+1;
if ~isempty(lambda_label)
    regularization.lambda = arginput{lambda_label};
    arginput(lambda_label-1:lambda_label) = [];
end

smooth_label = find(strcmpi(arginput,'smooth')==1)+1;
if ~isempty(smooth_label)
    smoothness.method = arginput{smooth_label};
    arginput(smooth_label-1:smooth_label) = [];
end

beta_label = find(strcmpi(arginput,'beta')==1)+1;
if ~isempty(beta_label)
    smoothness.beta = arginput{beta_label};
    arginput(beta_label-1:beta_label) = [];
end

conditionOn_label = find(strcmpi(arginput,'conditionOn')==1);
if ~isempty(conditionOn_label)
    conditionOn.variable = arginput{conditionOn_label+1};
    conditionOn.params   = arginput{conditionOn_label+2};
    
    index = cellfun(@(x) getfield(order, x), conditionOn.variable);
    indexEnd = cumsum(cell2mat(var.num_regress));
    indexStart = [1 indexEnd(1:end-1)+1];


    conditionOn.index = cell2mat(arrayfun(@(x) indexStart(x):indexEnd(x),...
        index,'UniformOutput', false));
    arginput(conditionOn_label:conditionOn_label+2) = [];
end

if isempty(arginput)
    sel_var = 1:length(var.name);
else
    for i = 1:length(arginput)
        eval(['sel_var(end+1) = order.',arginput{i},';']);
    end
end
smoothness.num_reg = cell2mat(var.num_regress(sel_var));

end


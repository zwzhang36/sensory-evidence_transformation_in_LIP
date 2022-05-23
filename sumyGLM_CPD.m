function F = sumyGLM_CPD(results_all, results_partial, p_criterion, tarVars)
%%%%%%%%%%%%%%%%%%%%%%%
%%% plot the CPDs
%%%%%%%%%%%%%%%%%%%%%%%
%%
signOnly = true;

index = zeros(size(results_all{1, 1}.label))==1;
for var = tarVars
    index(strcmp(results_all{1, 1}.label, var{:})) = true;
end

numVars  = numel(tarVars);

% significance test criterion
p_value  = p_criterion.p;
succ_len = p_criterion.successive_tw;
pvalue_str = repmat('1',1,succ_len);

% figure handle
F = {};
%%
coef_all = {};% coefficient for each variables
sign_time = []; % the first time showing a significance
for result = results_all'
    result = result{:};
    coef_all = [coef_all result.weight_mean];
    sign_time_curr = zeros(1,length(result.weight_mean));
    for ii = 1:length(result.weight_mean)
        % find the time becoming significance, succssive
        sign_ = num2str(result.pvalue{ii}'<p_value);
        sign_ = sign_(~isspace(sign_));
        
        if ~contains(sign_, pvalue_str)
            sign_time_curr(ii) = 100;
            pos_neg_curr(ii) = 0;
        else
            temp = strfind(sign_, pvalue_str);
            sign_time_curr(ii) = temp(1);
        end
    end
    sign_time = [sign_time; sign_time_curr];
end

% neurons with significant selevtivity
sign_neuron = arrayfun(@(x) find(sign_time(:,x)<100), 1:size(sign_time,2),...
    'UniformOutput', false);
sign_neuron = sign_neuron(index);
%% partial correlation of latencies among weight, consistent, evidence
% CPD = (SSE(reduced)-SSE(full))/SSE(reduced)

CPD = struct();
SSE_full = cellfun(@(x) x.SSE, results_all);

for varName = tarVars
    var = results_partial.(varName{:});
    SSE = cellfun(@(x) cellfun(@(y)  y.SSE, x, 'UniformOutput', false),...
        var, 'UniformOutput', false); % deviation
    SSE = cell2mat([SSE{:}]);
    CPD.(varName{:}) = mean((SSE - SSE_full)./SSE,2);
end

%% correlation between any two variable

% for i = 1:numVars
%     assert(strcmpi(varNames{i}, variable.name{i}))
% end

F{end+1} = figure('Name','correlation between CPDs','Units','normalized',...
    'Position',[0,0,0.35,0.4]);
for numVar1 = 1:numel(tarVars)
    for numVar2 = 1:numel(tarVars)
        if numVar1<=numVar2; continue; end
        
        if signOnly
            signNeu = intersect(sign_neuron{numVar1}, sign_neuron{numVar2});
        else
            signNeu = 1:nNeurons;
        end
        
        subplot(numVars-1,numVars-1,(numVar1-2)*(numVars-1)+numVar2);
        hold on;
        
        var1 = CPD.(tarVars{numVar1})(signNeu);
        var2 = CPD.(tarVars{numVar2})(signNeu);
        % orginal dots
        plot(var1, var2,'*');
        % fitted line
        lm = fitlm(var1, var2);
        [b,~,p] = gmregress(var1, var2);
        fplot(@(x) b(1) + b(2)*x, [min(var1),max(var1)]);
        text(max(var1), predict(lm, max(var1)),...
            sprintf('%.2e + %.2f*x\n p value:%.2e', b(1),...
            b(2), p(end)))
        % labels
        xlabel(tarVars{numVar1}); ylabel(tarVars{numVar2});
        
        % title(sprintf('pvalue: %0.5f', lm.Coefficients.pValue(end)));
    end
end

%
% subplot(numVars-1,numVars-1,3);
% axis off
title('correlation between CPDs of any two coefficients','FontSize', 14);

end



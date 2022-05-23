function F = sumyGLM_WCE(results, p_criterion)

%%%%%%%%%%%%%%%%%%%%%%%
%%% all the analysisi about single neuron selectivity in
%%% weight/consistency/evidence
%%% Author: Zhewei Zhang, 
%%% Date:   2020
%%%%%%%%%%%%%%%%%%%%%%%

%%

indMonkeyM = cellfun(@(x) strcmp(x.neuron(1), 'M'), results);
indMonkeyH = cellfun(@(x) strcmp(x.neuron(1), 'H'), results);

% variables that we are interested
vars = {'constant', 'weight', 'consis', 'tempevi', 'tempevi_color'};

ind = zeros(1, numel(results{1,1}.label));
for var = vars
    ind = ind | strcmp(results{1,1}.label, var{:});
end
  
for i = 1:numel(results)
    results{i}.pvalue = results{i}.pvalue(ind);
    results{i}.weight_std  = results{i}.weight_std(ind);
    results{i}.weight_mean = results{i}.weight_mean(ind);
end

ind(1) = [];
variable = results{1,1}.variable;
for field = fieldnames(variable)'
    variable.(field{:}) = variable.(field{:})(ind);
end
numVars  = numel(variable.name);

% time window of regression
whole_window = variable.ntfilt{1}*int2unit(variable.unit{1});
whole_regres = variable.num_regress{1};
regre_window = whole_window/whole_regres;
%%
% significance test criterion
p_value  = p_criterion.p;
succ_len = p_criterion.successive_tw;

pvalue_str = repmat('1',1,succ_len);
nNeurons   = numel(results);
%
signOnly  = false;

% figure handle
F = {};
%%
coef_all = {};% coefficient for each variables
sign_time = []; % the first time showing a significance
pos_neg = []; % the first time showing a significance
for result = results'
    result = result{:};
    coef_all = [coef_all result.weight_mean];
    sign_time_curr = zeros(1,length(result.weight_mean));
    pos_neg_curr = zeros(1,length(result.weight_mean));
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
            coef_sign = sign(result.weight_mean{ii}(sign_time_curr(ii):sign_time_curr(ii)+succ_len-1));
            if all(coef_sign>0)
                pos_neg_curr(ii) = 1;
            elseif all(coef_sign<0)
                pos_neg_curr(ii) = -1;
            else
                pos_neg_curr(ii) = 0;
            end
        end
    end
    sign_time = [sign_time; sign_time_curr];
    pos_neg = [pos_neg; pos_neg_curr];
end

%%
for index = 2: numel(vars)
    nloop = 0; eb = {};
    coef_all_temp = {coef_all(index, indMonkeyM), coef_all(index, indMonkeyH)};
    
    min_value = min(mean(cell2mat([coef_all_temp{:}]),2));
    Alpha = 0.01; minLen = 3; % it is confirmed by shuffle analysis
    F{end+1} = figure('Name', ['kernal_', vars{index}], 'Units', 'normalized',...
        'Position', [0,0,0.45,0.6]);  hold on
    colors = {[0.2 0.8 0.8], [1 1 0.4]};
    for coef = {coef_all_temp} %
        coef = coef{:};
        
        for i = 1:size(coef,2)
            nloop = nloop+1;
            coff_curr = cell2mat(coef{:,i});
            x = 1:size(coff_curr,1);
            y = mean(coff_curr,2);
            yerr = std(coff_curr,0,2)/sqrt(size(coff_curr,2));
            eb{end+1} = patch('XData', [x; flip(x)]', 'YData', [y-yerr; flip(y+yerr)]);
            set(eb{end}, 'FaceColor', colors{nloop}, 'EdgeColor', [0.2 0.2 0.2]);
            line(x, y,'Color', 0.5+(colors{nloop}-0.5)/2,'LineWidth',2)
            [h, ~] = ttest(coff_curr', 0, 'Alpha', Alpha);
            h = find(h);
            if any(h)
                % may multiple and seperate period with significance
                sega = find(diff(h)~=1);
                if ~isempty(sega)
                    sega = [sega(1) diff(sega) numel(h)-sega(end)];
                    significance = mat2cell(h,1,sega);
                else
                    significance = {h};
                end
                significance(cellfun(@numel, significance)<minLen) = [];
                % rectangle plot
                for ii = 1:numel(significance)
                    x_start = significance{ii}(1)-0.5;
                    y_hight = 0.003;
                    y_start = min_value-y_hight*i;
                    x_length = significance{ii}(end) + 0.5 - x_start;
                    rectangle('Position',[x_start y_start x_length y_hight],'FaceColor',...
                        eb{end}.FaceColor,'EdgeColor',[1 1 1])
                end
            end
        end
        
        ylabel('Gain');
        xlabel('time from the stimulus on(ms)');
        set(gca, 'FontSize', 12, ...
            'xtick',     0:whole_regres/3:whole_regres,...
            'xticklabel',num2str([0:whole_window/3:whole_window]') );
        legend([eb{:}], {'monkey M','monkey H'},'Location','southeast');
        title(['kernel_',vars{index}]);
    end
    
    plot([0 50], [0 0], 'k--', 'LineWidth', 1);
end
%%
% sort the neurons by the mean value of the selectivties to the evidence
% and plot the dynamics of the neurons' selectivity
pos = strcmpi('weight', results{1, 1}.label);
coeffTemp = cell2mat(coef_all(pos,:));

labels_M = []; labels_H = [];
if sum(indMonkeyM) ~= 0
    coeff_M = pca(coeffTemp(:, indMonkeyM));
    [~, labels_M] = sort(coeff_M(:,1));
end
if sum(indMonkeyH) ~= 0
    coeff_H = pca(coeffTemp(:, indMonkeyH));
    [~, labels_H] = sort(coeff_H(:,1));
end

if indMonkeyM(1)
%     labels = [labels_M; sum(indMonkeyM)+labels_H];
    labels = [sum(indMonkeyM)+labels_H; labels_M];
else 
    labels = [labels_H; sum(indMonkeyM)+labels_M];
end

F{end+1}= figure('Name','neural dynamics','Units','normalized',...
    'Position',[0,0,0.55,0.4]); 
for i = 2:1+numVars
    subplot(1,numVars,i-1)
    % coefficients for current variable
    coff_curr = cell2mat(coef_all(i,:));
    coff_sort = coff_curr(:,labels);
    
    % plot the heat map
    coff_sort = [coff_sort(1,:);coff_sort];
    x = 0:size(coff_sort,1)-1;
    y = 1:size(coff_sort,2);
    [xx,yy] = meshgrid(x,y);
    
    surf(xx,yy,coff_sort');
    shading flat
    view([0 0 1]);
    colormap(jet);
    %     caxis([-1,1]);
    colorbar;
    
    title([vars{i}]);
    xlabel 'time from the shape on(ms)'
    ylabel 'LIP neurons'
    axis([0 whole_regres 0.5 nNeurons]);
    set(gca,...
        'xtick',     0:whole_regres/3:whole_regres,...
        'xticklabel',num2str([0:whole_window/3:whole_window]') );
    colormapRange = caxis;
    caxis([-max(abs(colormapRange)), max(abs(colormapRange))])
end

%% plot the latency for each variables, only neurons with significance into account
F{end+1}= figure('Name','latencies','Units','normalized',...
    'Position',[0,0,0.45,0.6]); hold on;
% number of neurons showing significance
sign_num = sum(sign_time~=100,1);
% latency of significant neurons
latency_sign = mat2cell(sign_time(sign_time<100)',1, sign_num);

latency_sign_allvars = [];
title_ = [];
for i = 2:1+numVars
    % prepare latency data and their labels for box plot and multcompare
    labels = i*ones(size(latency_sign{i}));
    latency_sign_allvars = [latency_sign_allvars, [latency_sign{i};labels]];
    
    % sactter plot, latency for each neuron
    x = i-1.2+0.4*rand(sign_num(i),1);
    y = latency_sign{i}';
    title_ = [title_, sprintf('%s latency: %4.2f +- %4.2f\n', ...
            vars{i},...
            mean(y)*regre_window,...
            std(y)/sqrt(numel(y))*regre_window)];
    plot(x,y,'.k');
end
boxplot(latency_sign_allvars(1,:),latency_sign_allvars(2,:),'Widths',0.5)
set(gca,'xtick',1:numVars,'xticklabel',variable.name);
set(gca,'ytick',0:10:50,'yticklabel',{'0','300','600','900','1200','1500'});% TODO:
ylabel('latency (ms)');
[p,~,stats] = anova1(latency_sign_allvars(1,:),latency_sign_allvars(2,:),'off');
title_ = [title_, sprintf('p value: %4.2e', p)];
title(title_);
end


function  unit = int2unit(interval)
if ischar(interval)
    unit = interval;
    unit(end-1:end)=[];
    unit = str2double(unit);
    return
end
if isfloat(interval)
    unit = interval;
end
end
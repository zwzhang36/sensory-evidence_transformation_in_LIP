function F = sumyGLM_InOut(results, p_criterion, tw)

% GLM on the logLR(Tin) and logLR(Tout)
% reproduce the right figure in the slide 18 in the 20180521 proposal's ppt
% time window for single neuron summary
variable = results{1,1}.variable;
whole_window = variable.ntfilt{1}*int2unit(variable.unit{1});
whole_regres = variable.num_regress{1};

tw_start = ceil(tw(1)/whole_window*whole_regres)+1;
tw_end   = ceil(tw(2)/whole_window*whole_regres);
twInd    = tw_start:tw_end;

% significance test criterion
p_value  = p_criterion.p;
succ_len = p_criterion.successive_tw;

pvalue_str = repmat('1',1,succ_len);
nNeurons = numel(results);
%% find the time neurons with selectivity
coff = {};% coefficient for each variables
sign_time = {}; % the first time showing a significance
for result = results'
    result = result{:};
    coff = [coff result.weight_mean];
    sign_time_curr = cell(1,length(result.weight_mean));
    for ii = 1:length(result.weight_mean)
        % find the time becoming significance, succssive 
        if ii~=1
            sign = num2str(result.pvalue{ii}(twInd)'<p_value);
        else
            sign = num2str(result.pvalue{ii}'<p_value);
        end
        sign = sign(~isspace(sign));
        
        sign_time_curr{ii} = strfind(sign, pvalue_str);
        if isempty(sign_time_curr{ii})
            sign_time_curr{ii} = 100;
        else
            sign_time_curr{ii} = sign_time_curr{ii}(1);
        end
    end
    sign_time = [sign_time sign_time_curr'];
end
sign_time = cell2mat(sign_time)';

%% plot averaged kernal across all neurons
F = {};

vars = {'constant', 'tempevi_Tin', 'tempevi_Tout'};
F{end+1} = figure('Units', 'normalized', 'Position', [0,0,0.45,0.6]);  hold on
colors = {[0.2 0.2 0.8], [0.8 0.2 0.2]};

for index = 2: numel(vars)
    
    coff_curr = cell2mat(coff(index, :));
    
    x = 1:size(coff_curr,1);
    y = mean(coff_curr,2);
    yerr = std(coff_curr,0,2)/sqrt(size(coff_curr,2));
    eb = patch('XData', [x; flip(x)]', 'YData', [y-yerr; flip(y+yerr)]);
    set(eb, 'FaceColor', colors{index-1}, 'EdgeColor', [0.2 0.2 0.2]);
    line(x, y,'Color', 0.5+(colors{index-1}-0.5)/2,'LineWidth',2)
    
    ylabel('Gain');
    xlabel('time from the stimulus on(ms)');
    set(gca, 'FontSize', 12, ...
        'xtick',     0:whole_regres/3:whole_regres,...
        'xticklabel',num2str([0:whole_window/3:whole_window]') );
    title(['kernel_',vars{index}]);
    
    
    plot([0 60], [0 0], 'k--', 'LineWidth', 1);
end

%% scatter plot showing the coeffcient in 750ms/1000ms-1500/1750ms

Tin_pos  = find(strcmp(result.variable.name,'tempevi_Tin'))+1;
Tout_pos = find(strcmp(result.variable.name,'tempevi_Tout'))+1;

coef_Tin  = cell2mat(coff(Tin_pos,:));  coef_Tin  = coef_Tin(twInd, :);
coef_Tout = cell2mat(coff(Tout_pos,:)); coef_Tout = coef_Tout(twInd, :);

coef_Tin_mean  = mean(coef_Tin,1);
coef_Tout_mean = mean(coef_Tout,1);

[nodiff, sel_Tin ,sel_Tout] = deal(zeros(nNeurons,1));
for i = 1:nNeurons
    nodiff(i)  = ttest2(coef_Tin(:,i), -coef_Tout(:,i),'Alpha',0.01);
    sel_Tin(i)  = ttest(coef_Tin(:,i), 0,'Alpha',0.01);
    sel_Tout(i) = ttest(coef_Tout(:,i),0,'Alpha',0.01);
end


signNeuron = struct('both_nodiff', nodiff==0 & sel_Tin==1 & sel_Tout==1,...
                    'both_diff'  , nodiff==1 & sel_Tin==1 & sel_Tout==1,...
                    'InOnly'     , sel_Tin==1 & sel_Tout==0,...
                    'OutOnly'    , sel_Tin==0 & sel_Tout==1);
markers = struct('both_nodiff','^k', 'both_diff','^b', 'InOnly', '^r','OutOnly','^g');
labels  = struct('both_nodiff','logLR(Tin) == logLR(Tout)','InOnly', 'logLR(Tin)',...
                 'both_diff',  'logLR(Tin) != logLR(Tout)','OutOnly','logLR(Tout)');
%
F{end+1} = figure('Name','comparsion_TinTout'); 
axes1 = axes('Parent',F{end}); hold(axes1,'on');
plot(coef_Tin_mean, coef_Tout_mean,'^','MarkerSize',6,'DisplayName','no significance');
for type = fieldnames(signNeuron)'
    neuron = getfield(signNeuron, type{:});
    marker = getfield(markers,    type{:});  label = getfield(labels, type{:});
    plot(coef_Tin_mean(neuron),...
         coef_Tout_mean(neuron),...
         marker,'MarkerSize',6,'markerFaceColor',marker(end),...
         'DisplayName',label);
end
set(gca,'XAxisLocation','origin','YAxisLocation','origin');axis equal
legend(axes1,'show');

x = -0.2:0.01:0.3;y=-x; 
plot(x,y,'linewidth',2,'Parent',axes1)% diagonal line

x = -0.2:0.01:0.3;
fit = fitlm(coef_Tin_mean, coef_Tout_mean);
y= fit.Coefficients.Estimate(1,1)+fit.Coefficients.Estimate(2,1)*x;
plot(x,y,'linewidth',2,'Parent',axes1)% fitting line


%%  

F{end+1} = figure('Name','comparsion_TinTout_statis');
axes2 = axes('Parent',F{end}); hold(axes2,'on');
h1 = histogram(coef_Tin_mean+coef_Tout_mean,-0.2:0.05:0.45,  'FaceColor','r','DisplayName','logLR(Tin)',...
    'Parent',axes2,'Normalization','probability');
% h2 = histogram(-coef_Tout_mean,'FaceColor','g','DisplayName','logLR(Tout)',...
%     'Parent',axes2,'NumBins',10,'Normalization','probability');
% ymax = max([h1.Values, h2.Values]);
ymax = max([h1.Values]);
ylim([0, 1.5*ymax]);

% legend(axes2,'show');

% plot the mean value arrow 
pos = get(gca, 'Position');
color = {'r','g'};
position = {[ymax+0.15, ymax+0.07], [ymax+0.1, ymax+0.02]};

% Positins for the end of the Arrow in data units.
ith = 1;
for value = [coef_Tin_mean+coef_Tout_mean]'
    [xPos_start, xPos_end] = deal(mean(value));
    [yPos_start, yPos_end] = deal(position{ith}(1), position{ith}(2));
    % Create a textarrow annotation at the coordinates in data units
    % the textarrow coordinates are given [end_x, head_x], [end_y, head_y]
    x = [(xPos_start + abs(min(xlim)))/diff(xlim) * pos(3) + pos(1),...
        (xPos_end + abs(min(xlim)))/diff(xlim) * pos(3) + pos(1) ];
    y = [(yPos_start - min(ylim))/diff(ylim) * pos(4) + pos(2),...
        (yPos_end - min(ylim))/diff(ylim) * pos(4) + pos(2)];
    y(y>1)=1;y(y<0)=0;
    annotation(F{end},'textarrow',x,y,'String',[num2str(mean(value)),...
        '+',num2str(std(value))],'Color',color{ith});
    ith = ith + 1;
end
% significance test
[~,p] = kstest2(coef_Tin_mean,-coef_Tout_mean);
title(['ks test for distribution logLR(Tin)/(Tout) || p value:', num2str(p)]);
[~,p] = kstest2(coef_Tin_mean+coef_Tout_mean, 0);
title(['t test for distribution logLR(Tin)/(Tout) || p value:', num2str(p)]);
%% disp the number of each type of neurons
disp(['----------------------------------------------------']);
disp(['ks test for distribution logLR(Tin) and logLR(Tout)||p_value:', num2str(p)]);

disp(['logLR(Tin) == logLR(Tout):', num2str(sum(signNeuron.both_nodiff == 1))]);
disp(['logLR(Tin) != logLR(Tout):', num2str(sum(signNeuron.both_diff   == 1))]);
disp(['only logLR(Tin) != 0     :', num2str(sum(signNeuron.InOnly      == 1))]);
disp(['only logLR(Tout) != 0    :', num2str(sum(signNeuron.OutOnly     == 1))]);
disp(['----------------------------------------------------']);


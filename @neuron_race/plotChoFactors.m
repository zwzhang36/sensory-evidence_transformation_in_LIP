function [F, theta] = plotChoFactors(obj, varargin)
%% Logistic regression/ only for the shapes with the same number of shapes

monkey = obj.basic_info.monkey;
race = obj.race;
tarpos = race.tarpos;
choice = race.choice_color;
condition = cell2mat(race.condition');

if contains(monkey, 'M')
    reward = 1-race.rewardfeedback;
elseif contains(monkey, 'H')
    reward = race.outcome;
end

if ~isempty(varargin) && strcmpi(varargin{1}, 'pre-assigned')
    weightList = [0.9 0.54 0.3 -0.3 -0.54 -0.9 0.9 0.54 0.3 -0.3 -0.54 -0.9];
else
    weightList = [obj.subweight(1:6), -obj.subweight(7:end)];
end
%%
Red    = condition < 7;
Gre    = condition >=7;
Weight = weightList(condition);

Weight_Red = Weight.*Red;
Weight_Gre = Weight.*Gre;

tarpos = 2*tarpos-1; % 1:red left; -1;red right
reward = 2*reward-1; % monkey M: 1/0:long/short reward; H: 1/0:long,short/no reward

red_long = choice.*reward;
gre_long = (~choice).*reward;

%% regression

X = Weight_Red - Weight_Gre;
Y = choice';

% logistic regression
[theta, ~, stats] = glmfit(X,Y,'binomial','link','logit','constant','on');

%% plot
F = figure('Name','temporal effect','visible','on','Units','normalized',...
    'Position',[0,0,0.4,0.5]);hold on;

theta = theta(2:end);
bar(theta, 0.7);
errorbar(1:numel(theta), theta, stats.se(2:end), '.','MarkerSize',10)

title 'Logistic regression'
set(gca,'FontSize',14);
ylabel 'Coefficients'

ylim([0, 1.25]); xlim([0.5 6.5])

set(gca, 'xtick',1:6,'xticklabel',{'1st','2nd','3th','4th','5th','6th'},...
    'fontsize', 12);

end